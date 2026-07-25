#!/usr/bin/env python3
"""Test a named set of Y markers against a BAM at read level.

Motivation
----------
`y_haplo_from_markers.py` selects markers by their ISOGG label, which is fine for
the ISOGG-covered backbone but blind to everything ISOGG never ranked -- markers
carrying ``isogg_haplogroup=unknown`` / ``not listed``. FTDNA's Block Tree and
YFull both name such markers routinely, so a block of SNP names copied from
either tool could not previously be tested at all.

This script takes marker *names*, resolves them through the index built by
``annotate/build_marker_index.sh``, and reports ancestral/derived read counts
straight from the BAM. It applies no label filtering.

Example
-------
    annotate/y_markers_pileup.py \
        --bam sample.bam --ref hg38.fa \
        --markers Z6128,S17295,FT156872 \
        --out results/block_test.tsv

    annotate/y_markers_pileup.py --bam sample.bam --ref hg38.fa \
        --marker-file ftdna_block.txt --label G-PF3239
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter

INDEL = re.compile(r"[+-](\d+)")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
DEAMINATION = {("C", "T"), ("G", "A")}


def parse_bases(s: str, ref: str) -> Counter:
    """Count called bases in an mpileup base string."""
    out: Counter = Counter()
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "^":  # read start; next char encodes MAPQ
            i += 2
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            m = INDEL.match(s, i)
            if not m:  # malformed; skip the marker char
                i += 1
                continue
            i = m.end() + int(m.group(1))
            continue
        if ch in ".,":
            out[ref.upper()] += 1
        elif ch.upper() in "ACGTN":
            out[ch.upper()] += 1
        elif ch in "*#":
            out["del"] += 1
        i += 1
    return out


def mutation_class(anc: str, der: str) -> str:
    if len(anc) != 1 or len(der) != 1:
        return "indel/other"
    if (anc, der) in DEAMINATION:
        return "transition(deamination-prone)"
    if (anc, der) in TRANSITIONS:
        return "transition"
    return "transversion"


def site_call(dp: int, a: int, d: int, mclass: str) -> str:
    """Pre-registered call rules. Must stay identical to y_sites_pileup.py.

    Read counts alone are not a call. A single read is one molecule, and at aDNA
    depths a single deamination-prone C>T/G>A read is the expected artifact
    rather than evidence -- documented in this repo for CGG017682, whose entire
    apparent G-panel signal was single-read transitions.

    The asymmetry that matters most: one ancestral read is `low_power`, never
    `ancestral`. Calling it ancestral converts "we could not tell" into "we
    tested and it was negative", which is the exact confusion the accompanying
    pre-registration exists to prevent.
    """
    if dp == 0:
        return "no_coverage"
    if d >= 2 and a == 0:
        return "DERIVED"
    if d == 1 and a == 0:
        if mclass == "transversion":
            return "DERIVED_1read_transversion"
        if mclass == "transition(deamination-prone)":
            return "nocall_damage_prone_1read"
        return "nocall_1read_transition"
    if a >= 2 and d == 0:
        return "ancestral"
    if a == 1 and d == 0:
        return "low_power_1read_ancestral"
    if a and d:
        return "mixed"
    return "other_allele"


def load_markers(index_path: str, wanted: set[str]) -> dict[str, dict]:
    found: dict[str, dict] = {}
    opener = gzip.open if index_path.endswith(".gz") else open
    with opener(index_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(header)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            name = f[col["name"]]
            if name in wanted:
                found[name] = {
                    "chrom": f[col["chrom"]],
                    "pos": int(f[col["pos"]]),
                    "anc": f[col["anc"]],
                    "der": f[col["der"]],
                    "isogg": f[col["isogg"]],
                    "yfull_node": f[col["yfull_node"]],
                }
    return found


def run_pileup(bam, ref, markers, min_mq, min_bq, flank, max_depth, workdir):
    """Extract only the reads overlapping the marker positions, then pile up.

    Piling up the whole contig is orders of magnitude slower -- samtools walks
    every base of chrY regardless of ``-l`` -- so slice the reads out first.
    """
    bed = os.path.join(workdir, "targets.bed")
    posf = os.path.join(workdir, "targets.pos")
    sub = os.path.join(workdir, "targets.bam")

    rows = sorted({(m["chrom"], m["pos"]) for m in markers.values()})
    with open(bed, "w") as fh:
        for chrom, pos in rows:
            fh.write(f"{chrom}\t{max(0, pos - 1 - flank)}\t{pos + flank}\n")
    with open(posf, "w") as fh:
        for chrom, pos in rows:
            fh.write(f"{chrom}\t{pos}\n")

    subprocess.run(
        ["samtools", "view", "-b", "-M", "-L", bed, bam, "-o", sub],
        check=True,
    )
    subprocess.run(["samtools", "index", sub], check=True)

    proc = subprocess.run(
        [
            "samtools", "mpileup", "-f", ref, "-l", posf,
            "-q", str(min_mq), "-Q", str(min_bq),
            "-d", str(max_depth), "--no-BAQ", sub,
        ],
        check=True, capture_output=True, text=True,
    )

    pile = {}
    for ln in proc.stdout.splitlines():
        c = ln.split("\t")
        if len(c) < 5:
            continue
        pile[(c[0], int(c[1]))] = (c[2], int(c[3]), parse_bases(c[4], c[2]))
    return pile


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--index", default=None,
                    help="marker index (default: resources/marker_index.tsv.gz)")
    ap.add_argument("--markers", default=None, help="comma-separated marker names")
    ap.add_argument("--marker-file", default=None,
                    help="file of marker names, one per line (blank lines and # ignored)")
    ap.add_argument("--label", default="", help="block/node label recorded in the output")
    ap.add_argument("--min-mq", type=int, default=25)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--flank", type=int, default=300,
                    help="bp of read-extraction padding around each marker")
    ap.add_argument("--max-depth", type=int, default=1000)
    ap.add_argument("--out", default="-")
    args = ap.parse_args()

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    index_path = args.index or os.path.join(repo_root, "resources", "marker_index.tsv.gz")
    if not os.path.exists(index_path):
        sys.exit(f"marker index not found: {index_path}\n"
                 f"build it with: annotate/build_marker_index.sh")

    names: list[str] = []
    if args.markers:
        names += [x.strip() for x in args.markers.split(",") if x.strip()]
    if args.marker_file:
        with open(args.marker_file) as fh:
            for ln in fh:
                ln = ln.split("#", 1)[0].strip()
                if ln:
                    names.append(ln)
    if not names:
        sys.exit("no markers given (use --markers or --marker-file)")

    # de-duplicate, preserve input order
    seen: set[str] = set()
    ordered = [n for n in names if not (n in seen or seen.add(n))]

    markers = load_markers(index_path, set(ordered))
    missing = [n for n in ordered if n not in markers]
    if missing:
        print(f"# not in catalogue ({len(missing)}): {','.join(missing)}", file=sys.stderr)
    if not markers:
        sys.exit("none of the requested markers are in the catalogue")

    with tempfile.TemporaryDirectory(prefix="ymarkers.") as workdir:
        pile = run_pileup(args.bam, args.ref, markers, args.min_mq, args.min_bq,
                          args.flank, args.max_depth, workdir)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    cols = ["block", "marker", "chrom", "pos", "anc>der", "mutation_class",
            "depth", "anc_reads", "der_reads", "other_reads", "call",
            "isogg", "yfull_node"]
    print("\t".join(cols), file=out)

    tally: Counter = Counter()
    for name in ordered:
        m = markers.get(name)
        if m is None:
            print("\t".join([args.label, name, ".", ".", ".", ".", "0", "0", "0", "0",
                             "not_in_catalogue", ".", "."]), file=out)
            tally["not_in_catalogue"] += 1
            continue
        _ref, dp, cnt = pile.get((m["chrom"], m["pos"]), (".", 0, Counter()))
        a = cnt.get(m["anc"], 0)
        d = cnt.get(m["der"], 0)
        other = sum(v for k, v in cnt.items() if k not in (m["anc"], m["der"]))
        mclass = mutation_class(m["anc"], m["der"])
        call = site_call(dp, a, d, mclass)
        tally[call] += 1
        print("\t".join([
            args.label, name, m["chrom"], str(m["pos"]),
            f"{m['anc']}>{m['der']}", mclass,
            str(dp), str(a), str(d), str(other), call,
            m["isogg"], m["yfull_node"],
        ]), file=out)

    if out is not sys.stdout:
        out.close()

    print("# " + ", ".join(f"{k}={v}" for k, v in sorted(tally.items())), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
