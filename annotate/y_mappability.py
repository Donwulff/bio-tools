#!/usr/bin/env python3
"""Measure, per marker, whether our own pipeline can recover a read at all.

Why this exists
---------------
Every call this repository makes rests on reads that ``bwa aln`` chose to place
at a marker. Nothing in the pileup can see the reads that were placed somewhere
*else* -- ``ylib.mapq_audit()`` reports MAPQ at the site, so it sees ambiguous
reads that stayed and is structurally blind to reads that left. A site can
therefore look clean at ``0% MQ0`` while quietly losing most of its molecules to
an alt contig, a paralog, or a different build's version of the same region.

This tool measures the loss directly. For each marker it cuts the source
reference into every read of length L that overlaps that position -- L reads at
1 bp step -- and maps them back with the pipeline's own aligner and options. A
read that fails to return is a pure mappability failure: the reads are exact
copies of the reference, so there is no sequencing error, no damage and no
biology to disentangle. That is deliberate. Simulating errors would only blur
the measurement; damage belongs in a separate pass if its marginal effect is
wanted.

No liftover is required, and that is the point of the design. The tile carries
its own truth: L reads spanning one position should land on L consecutive
positions in the target, so the *modal* implied position is the marker's
location in that build, discovered rather than asserted. This makes the tool
work against CHM13 or any other assembly for which no chain file is held, and it
yields the cross-build coordinate as a by-product.

What the numbers mean
---------------------
``frac_recovered`` is the fraction of possible reads that come back to the right
locus at MAPQ >= 25 -- the same threshold the pileup applies. It is a power
figure: at 0.60, a site with nominal 10x coverage is really being interrogated
at 6x, and a marker at 0.00 is uncallable at that read length no matter how deep
the library. ``ref_matches_anc`` is a separate, cheaper trap: if a build's
reference carries the *derived* allele, ref and alt are swapped there and any
naive call against it is inverted. CHM13v2's chrY is a different individual's
chromosome from GRCh38's, so this is not hypothetical.

Usage
-----
    annotate/y_mappability.py \\
        --markers markers/L166_defining.txt markers/backbone_control.txt \\
        --source mapping/index/hg38p14DH3630O.fa \\
        --target working=mapping/index/hg38p14DH3630O.fa \\
        --target noalt=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \\
        --target chm13=mapping/index/chm13v2.0_maskedY_rCRSDH3630O.fa \\
        --target hs37d5=/mnt/GenomicData/hs37d5/bwa/hs37d5.fa \\
        --read-lengths 35,45,60,90 \\
        --out results/mappability/y_marker_mappability.tsv

The source reference supplies the truth coordinates and must be the build our
marker index is in (GRCh38). Targets are ``label=path`` and are processed one at
a time: bwa aln holds the whole index resident, so parallelism belongs in
--threads, not in running several references at once.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import os
import re
import subprocess
import sys
import tempfile

DEFAULT_INDEX = "resources/marker_index.tsv.gz"
DEFAULT_ALN_OPTS = "-n 0.01 -k 2 -l 1024"   # must match mapping/map_se_adna.sh
MIN_MQ = 25                                  # must match the pileup threshold
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_marker_index(path: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with _open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= len(hdr):
                r = dict(zip(hdr, f))
                out[r["name"]] = r
    return out


def load_marker_names(paths: list[str]) -> list[str]:
    names, seen = [], set()
    for p in paths:
        with open(p) as fh:
            for ln in fh:
                ln = ln.split("#", 1)[0].strip()
                if ln and ln not in seen:
                    seen.add(ln)
                    names.append(ln)
    return names


def faidx(ref: str, region: str) -> str:
    """Fetch a reference slice as a bare uppercase sequence."""
    p = subprocess.run(["samtools", "faidx", ref, region],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return ""
    return "".join(p.stdout.split("\n")[1:]).upper()


def ref_base(ref: str, chrom: str, pos: int) -> str:
    s = faidx(ref, f"{chrom}:{pos}-{pos}")
    return s if s else "?"


def substitute(seq: str, length: int, base: str) -> str:
    """Put `base` at the marker position (the centre of a 2L-1 window).

    Tiling straight from the reference gives every read the ANCESTRAL allele, so
    a recovery figure obtained that way measures ancestral-read mappability only.
    A read carrying the derived allele scores one mismatch lower against this
    locus and may score better against a paralogous one -- which is a different
    question and, for a marker two curators have excluded or down-rated, the more
    important one.
    """
    return seq[:length - 1] + base + seq[length:]


def build_tile(seq: str, marker: str, length: int) -> list[tuple[str, str]]:
    """Every read of `length` overlapping the centre of `seq`.

    `seq` is the 2L-1 window centred on the marker, so read i (0-based) places
    the marker at offset L-1-i from its start. Reads containing N are dropped --
    an N is not a mappability failure, it is an absent reference base, and
    counting it as failure would understate recovery.
    """
    reads = []
    for i in range(length):
        sub = seq[i:i + length]
        if len(sub) != length or "N" in sub:
            continue
        offset = length - 1 - i
        reads.append((f"{marker}|{length}|{offset}", sub))
    return reads


def norm_contig(c: str) -> str:
    """Compare contig names across builds that disagree about the 'chr' prefix.

    GRCh38/CHM13 write 'chrY'; hs37d5 and other GRCh37 builds write 'Y'. Without
    this every read mapped to a GRCh37 target counted as off-target and recovery
    read 0.000 for the whole reference -- a tool bug that looks exactly like the
    catastrophic finding it is not.
    """
    return c[3:] if c.startswith("chr") else c


def cigar_is_exact(cig: str) -> bool:
    """True for a clean ungapped alignment: a single M/=/X run, no clipping."""
    ops = CIGAR_RE.findall(cig)
    return len(ops) == 1 and ops[0][1] in "M=X"


def map_reads(fq: str, ref: str, threads: int, aln_opts: str,
              workdir: str) -> list[list[str]]:
    """bwa aln + samse with the pipeline's own options; returns SAM fields."""
    sai = os.path.join(workdir, "tile.sai")
    with open(sai, "wb") as out, open(os.devnull, "wb") as null:
        subprocess.run(["bwa", "aln", "-t", str(threads), *aln_opts.split(),
                        ref, fq], stdout=out, stderr=null, check=True)
    p = subprocess.run(["bwa", "samse", ref, sai, fq],
                       capture_output=True, text=True)
    os.unlink(sai)
    rows = []
    for ln in p.stdout.split("\n"):
        if ln and not ln.startswith("@"):
            rows.append(ln.split("\t"))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markers", required=True, nargs="+")
    ap.add_argument("--index", default=DEFAULT_INDEX)
    ap.add_argument("--source", required=True,
                    help="reference the marker index is in (GRCh38); supplies truth")
    ap.add_argument("--target", required=True, action="append",
                    metavar="LABEL=PATH", help="repeatable; bwa-indexed reference")
    ap.add_argument("--read-lengths", default="45",
                    help="comma-separated, e.g. 35,45,60,90")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--aln-opts", default=DEFAULT_ALN_OPTS)
    ap.add_argument("--min-mq", type=int, default=MIN_MQ)
    ap.add_argument("--allele", choices=("anc", "der"), default="anc",
                    help="which allele the tiled reads carry at the marker "
                         "position (default anc, i.e. straight from the "
                         "reference). 'der' tests whether a derived read maps "
                         "back to the same locus -- allele asymmetry is a "
                         "reliability failure a reference-only tile cannot see.")
    ap.add_argument("--out", default="-")
    a = ap.parse_args()

    lengths = [int(x) for x in a.read_lengths.split(",") if x.strip()]
    targets = []
    for t in a.target:
        if "=" not in t:
            print(f"--target needs LABEL=PATH, got {t!r}", file=sys.stderr)
            return 1
        label, path = t.split("=", 1)
        if not os.path.exists(path + ".bwt"):
            print(f"no bwa index beside {path}", file=sys.stderr)
            return 1
        targets.append((label, path))

    idx = load_marker_index(a.index)
    names = load_marker_names(a.markers)
    markers = []
    for n in names:
        r = idx.get(n)
        if r is None:
            print(f"SKIP {n}: not in marker index", file=sys.stderr)
            continue
        markers.append((n, r["chrom"], int(r["pos"]), r["anc"], r["der"]))
    if not markers:
        print("no markers resolved", file=sys.stderr)
        return 1

    out = sys.stdout if a.out == "-" else open(a.out, "w")
    cols = ["marker", "chrom_src", "pos_src", "anc", "der", "allele", "read_len", "target",
            "n_reads", "n_mapped", "n_exact_correct_contig", "n_mq_ge_min",
            "frac_recovered", "n_unique", "frac_recovered_unique",
            "n_mq0", "n_offtarget", "modal_pos",
            "modal_support", "target_ref_base", "ref_matches_anc",
            "top_offtarget"]
    print("\t".join(cols), file=out)

    work = tempfile.mkdtemp(prefix="ymap.")
    try:
        # One FASTQ for every (marker, read length): bwa aln holds a multi-GB
        # index resident and loading it dominates a job this small, so each
        # reference is loaded exactly once. The read name carries its own length,
        # so the results separate again on the way out.
        tiles, per_marker = [], {}
        for length in lengths:
            for name, chrom, pos, anc, der in markers:
                seq = faidx(a.source, f"{chrom}:{pos-length+1}-{pos+length-1}")
                if len(seq) != 2 * length - 1:
                    print(f"SKIP {name} at L={length}: short/absent window",
                          file=sys.stderr)
                    continue
                if a.allele == "der":
                    seq = substitute(seq, length, der)
                t = build_tile(seq, name, length)
                per_marker[(name, length)] = len(t)
                tiles.extend(t)
        if not tiles:
            print("no tiles built", file=sys.stderr)
            return 1

        fq = os.path.join(work, "tile.fq")
        with open(fq, "w") as fh:
            for rid, seq in tiles:
                fh.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")

        src_chrom_of = {m[0]: m[1] for m in markers}
        for label, path in targets:
            rows = map_reads(fq, path, a.threads, a.aln_opts, work)
            agg = collections.defaultdict(lambda: {
                "mapped": 0, "exact": 0, "mq": 0, "uniq": 0, "mq0": 0, "off": 0,
                "pos": collections.Counter(), "offc": collections.Counter()})
            for f in rows:
                if len(f) < 6:
                    continue
                qname, flag, rname, spos, mapq, cig = (
                    f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5])
                marker, lstr, off = qname.split("|")
                length, off = int(lstr), int(off)
                s_ = agg[(marker, length)]
                if flag & 4 or rname == "*":
                    continue
                s_["mapped"] += 1
                if mapq == 0:
                    s_["mq0"] += 1
                # The source contig is the truth; a read landing anywhere else
                # has been pulled away, which is exactly the loss a chrY-only
                # pileup can never see.
                if norm_contig(rname) != norm_contig(src_chrom_of[marker]):
                    s_["off"] += 1
                    s_["offc"][rname] += 1
                    continue
                if not cigar_is_exact(cig):
                    continue
                s_["exact"] += 1
                eff = (length - 1 - off) if (flag & 16) else off
                s_["pos"][spos + eff] += 1
                if mapq >= a.min_mq:
                    s_["mq"] += 1
                # Reported, not applied. `n_mq_ge_min` above is unchanged and
                # remains the criterion every call in this repository uses; this
                # counts what a uniqueness filter *would* have kept, so the two
                # can be compared on tiles whose truth is known by construction.
                # Adopting it is registered separately in
                # PREREG_uniqueness_filter.md and is not done here.
                xt = next((x[5:] for x in f[11:] if x.startswith("XT:A:")), None)
                x0 = next((int(x[5:]) for x in f[11:] if x.startswith("X0:i:")),
                          None)
                if xt == "U" and x0 == 1:
                    s_["uniq"] += 1

            for length in lengths:
                for name, chrom, pos, anc, der in markers:
                    n = per_marker.get((name, length))
                    if not n:
                        continue
                    s_ = agg[(name, length)]
                    modal_pos, modal_sup = (s_["pos"].most_common(1) or
                                            [(0, 0)])[0]
                    rb = "?"
                    if modal_sup:
                        # Contig naming differs across builds ("chrY" vs "Y"),
                        # so the base lookup tries the plausible spellings
                        # rather than reporting a missing base as a mismatch.
                        for cand in (chrom, chrom.replace("chr", ""), "chrY", "Y"):
                            rb = ref_base(path, cand, modal_pos)
                            if rb != "?":
                                break
                    top_off = ",".join(f"{k}:{v}" for k, v in
                                       s_["offc"].most_common(3)) or "-"
                    print("\t".join(str(x) for x in [
                        name, chrom, pos, anc, der, a.allele, length, label,
                        n, s_["mapped"], s_["exact"], s_["mq"],
                        f"{s_['mq'] / n:.3f}",
                        s_["uniq"], f"{s_['uniq'] / n:.3f}",
                        s_["mq0"], s_["off"],
                        modal_pos or "NA", modal_sup, rb,
                        "yes" if rb == anc else ("no" if rb != "?" else "NA"),
                        top_off]), file=out)
            out.flush()
            print(f"  {label}: {len(rows)} tile reads mapped", file=sys.stderr)
    finally:
        for f in os.listdir(work):
            os.unlink(os.path.join(work, f))
        os.rmdir(work)
        if out is not sys.stdout:
            out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
