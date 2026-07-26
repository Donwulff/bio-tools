#!/usr/bin/env python3
"""Lift this project's hg38 Y markers down to GRCh37, for pileup on a b37 deposit.

Why this exists
---------------
Every marker set here is in hg38, but some public ancient deposits are only
distributed as GRCh37/hs37d5 BAMs and are too shallow to be worth re-mapping.
Genotyping them means moving 20-odd coordinates down a build rather than moving
27 million reads up one.

Doing that by hand is where silent errors live, and there are three of them:

* **The coordinate.** chrY was substantially revised between builds, so an
  offset that holds for one marker does not hold for its neighbour. The chain
  file is read directly -- no ``liftOver`` binary -- and a marker landing in no
  chain block, or in an inverted one, is reported ``unmapped`` rather than
  given a plausible-looking wrong coordinate.
* **The contig name.** hg38 analysis sets call it ``chrY``; hs37d5 calls it
  ``Y``. Emitting the wrong one produces an empty pileup that looks exactly
  like no coverage, which is the single most dangerous failure mode in a test
  whose expected result is already "no coverage".
* **The reference base.** A lifted position is only trustworthy if the target
  build actually carries the ancestral allele there. This is checked against
  the target FASTA and reported per marker as ``ref_ok`` / ``REF_MISMATCH``.
  A mismatch means the lift landed somewhere wrong, or the marker's polarity
  differs between builds; either way the position must not be genotyped.

The output doubles as the ``--sites`` input for ``y_sites_pileup.py``, so the
coordinates that were checked are the coordinates that get used.

Usage
-----
    annotate/y_lift_markers.py \
        --markers markers/yfull_L166_defining.txt \
        --index resources/marker_index.tsv.gz \
        --chain /mnt/GenomicData/OpenSNP/puller/hg38ToHg19.over.chain \
        --target-ref /mnt/GenomicData/hs37d5/hs37d5.fa \
        --target-chrom Y \
        --sites-out sites_b37.tsv --report-out lift_report.tsv
"""

import argparse
import gzip
import subprocess
import sys

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from panel_membership import load_chain, lift  # noqa: E402


def load_index(path: str) -> dict[str, dict]:
    op = gzip.open if path.endswith(".gz") else open
    out: dict[str, dict] = {}
    with op(path, "rt") as fh:
        col = {name: i for i, name in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            out[f[col["name"]]] = {
                "chrom": f[col["chrom"]],
                "pos": int(f[col["pos"]]),
                "anc": f[col["anc"]],
                "der": f[col["der"]],
            }
    return out


def load_names(paths: list[str]) -> list[str]:
    names: list[str] = []
    for p in paths:
        with open(p) as fh:
            for ln in fh:
                ln = ln.split("#")[0].strip()
                if ln:
                    names.append(ln)
    return names


def ref_base(ref: str, chrom: str, pos1: int) -> str:
    """One base from the target FASTA, 1-based, uppercased."""
    out = subprocess.run(
        ["samtools", "faidx", ref, f"{chrom}:{pos1}-{pos1}"],
        capture_output=True, text=True, check=True,
    ).stdout.splitlines()
    return "".join(out[1:]).strip().upper() if len(out) > 1 else ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markers", nargs="+", required=True)
    ap.add_argument("--index", default="resources/marker_index.tsv.gz")
    ap.add_argument("--chain", required=True, help="hg38 -> hg19 chain file")
    ap.add_argument("--source-chrom", default="chrY")
    ap.add_argument("--target-chrom", default="Y",
                    help="contig name in the TARGET bam/reference (hs37d5 uses 'Y')")
    ap.add_argument("--target-ref", required=True)
    ap.add_argument("--sites-out", required=True)
    ap.add_argument("--report-out", required=True)
    args = ap.parse_args()

    index = load_index(args.index)
    names = load_names(args.markers)
    blocks = load_chain(args.chain, args.source_chrom)
    if not blocks:
        sys.exit(f"no same-strand {args.source_chrom} blocks in {args.chain}")

    rows, n_ok = [], 0
    for name in names:
        rec = index.get(name)
        if rec is None:
            rows.append((name, "NA", "NA", "NA", "NA", "NA", "not_in_catalogue"))
            continue
        # chain blocks are 0-based half-open; marker positions are 1-based.
        tgt0 = lift(rec["pos"] - 1, blocks)
        if tgt0 is None:
            rows.append((name, rec["pos"], "NA", rec["anc"], rec["der"], "NA", "unmapped"))
            continue
        tgt1 = tgt0 + 1
        rb = ref_base(args.target_ref, args.target_chrom, tgt1)
        status = "ref_ok" if rb == rec["anc"].upper() else f"REF_MISMATCH(found_{rb or 'none'})"
        if status == "ref_ok":
            n_ok += 1
        rows.append((name, rec["pos"], tgt1, rec["anc"], rec["der"], rb, status))

    with open(args.report_out, "w") as fh:
        fh.write("marker\tpos_hg38\tpos_target\tanc\tder\ttarget_ref_base\tstatus\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # Only ref_ok markers reach the pileup. A REF_MISMATCH is a failed lift, and
    # genotyping it would produce a confident call at the wrong place.
    with open(args.sites_out, "w") as fh:
        for name, _p38, ptgt, anc, der, _rb, status in rows:
            if status == "ref_ok":
                fh.write(f"{args.target_chrom}\t{ptgt}\t{anc}\t{der}\t{name}\n")

    print(f"# {len(names)} markers: {n_ok} lifted and ref-checked, "
          f"{len(names) - n_ok} excluded", file=sys.stderr)
    for r in rows:
        if r[6] != "ref_ok":
            print(f"#   excluded {r[0]}: {r[6]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
