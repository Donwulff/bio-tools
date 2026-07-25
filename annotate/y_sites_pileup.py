#!/usr/bin/env python3
"""Genotype arbitrary chrY positions with a mandatory MAPQ audit.

Companion to annotate/y_markers_pileup.py, which resolves *named* catalogue
markers. This one takes bare coordinates, which is what uncatalogued ("novel")
candidate sites necessarily are.

Allele counting uses the same thresholds as the marker script (-q 25 -Q 20
-d 1000 --no-BAQ). The MAPQ columns are computed separately from the raw
alignments, because depth alone cannot distinguish a real site from a
collapsed-repeat pileup: sites with unremarkable DP have been rejected at
39-88% MQ0.

Site file: TSV, no header, columns  chrom  pos  anc  der  [class]
Output:    TSV to stdout with one row per site.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter

INDEL = re.compile(r"[+-](\d+)")
DEAMINATION = {("C", "T"), ("G", "A")}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def parse_bases(s: str, ref: str) -> Counter:
    """Count called bases in an mpileup base string, tracking strand."""
    out: Counter = Counter()
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "^":          # read start, next char encodes MAPQ
            i += 2
            continue
        if ch == "$":          # read end
            i += 1
            continue
        if ch in "+-":         # indel, skip its payload
            m = INDEL.match(s, i)
            if m:
                i = m.end() + int(m.group(1))
                continue
            i += 1
            continue
        if ch in ".,":         # matches reference
            out[(ref.upper(), "+" if ch == "." else "-")] += 1
        elif ch.upper() in "ACGTN":
            out[(ch.upper(), "+" if ch.isupper() else "-")] += 1
        i += 1
    return out


def mapq_audit(bam: str, chrom: str, pos: int) -> tuple[int, int, int]:
    """Return (n_reads, n_mq0, n_mq60) for reads overlapping pos, unfiltered."""
    p = subprocess.run(["samtools", "view", bam, f"{chrom}:{pos}-{pos}"],
                       capture_output=True, text=True, check=True)
    n = n0 = n60 = 0
    for ln in p.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 5:
            continue
        flag = int(f[1])
        if flag & 0x900:       # skip secondary / supplementary
            continue
        mq = int(f[4])
        n += 1
        if mq == 0:
            n0 += 1
        if mq >= 60:
            n60 += 1
    return n, n0, n60


def pileup(bam: str, ref: str, chrom: str, pos: int,
           min_mq: int, min_bq: int, max_depth: int):
    p = subprocess.run(
        ["samtools", "mpileup", "-f", ref, "-r", f"{chrom}:{pos}-{pos}",
         "-q", str(min_mq), "-Q", str(min_bq), "-d", str(max_depth),
         "--no-BAQ", bam],
        capture_output=True, text=True, check=True)
    for ln in p.stdout.splitlines():
        c = ln.split("\t")
        if len(c) >= 5 and int(c[1]) == pos:
            return c[2].upper(), int(c[3]), parse_bases(c[4], c[2])
    return None, 0, Counter()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sites", required=True, help="TSV: chrom pos anc der [class]")
    ap.add_argument("--sample", default="")
    ap.add_argument("--min-mq", type=int, default=25)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--max-depth", type=int, default=1000)
    a = ap.parse_args()

    cols = ["sample", "chrom", "pos", "anc", "der", "mut_class", "dp", "n_anc",
            "n_der", "n_other", "fwd", "rev", "n_reads", "pct_mq0", "n_mq60", "call"]
    print("\t".join(cols))

    for ln in open(a.sites):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        chrom, pos, anc, der = f[0], int(f[1]), f[2].upper(), f[3].upper()
        mclass = f[4] if len(f) > 4 else (
            "transition" if (anc, der) in TRANSITIONS else "TRANSVERSION")

        _, dp, counts = pileup(a.bam, a.ref, chrom, pos,
                               a.min_mq, a.min_bq, a.max_depth)
        n_reads, n_mq0, n_mq60 = mapq_audit(a.bam, chrom, pos)

        n_anc = counts[(anc, "+")] + counts[(anc, "-")]
        n_der = counts[(der, "+")] + counts[(der, "-")]
        n_other = dp - n_anc - n_der
        fwd = sum(v for (b, s), v in counts.items() if s == "+")
        rev = sum(v for (b, s), v in counts.items() if s == "-")

        # Pre-registered decision rules. A single transition read is never
        # sufficient; single-read C>T/G>A is treated as deamination, not evidence.
        damage_prone = (anc, der) in DEAMINATION
        if dp == 0:
            call = "nocall_nocoverage"
        elif n_der >= 2 and n_anc == 0:
            call = "DERIVED"
        elif n_der == 1 and n_anc == 0 and not damage_prone and mclass == "TRANSVERSION":
            call = "DERIVED_1read_transversion"
        elif n_der == 1 and n_anc == 0 and damage_prone:
            call = "nocall_damage_prone_1read"
        elif n_anc >= 2 and n_der == 0:
            call = "ancestral"
        elif n_anc == 1 and n_der == 0:
            call = "low_power_1read_ancestral"
        elif n_anc and n_der:
            call = "MIXED_check_paralogy"
        else:
            call = "nocall"

        pct_mq0 = f"{100.0*n_mq0/n_reads:.0f}%" if n_reads else "NA"
        print("\t".join(str(x) for x in [
            a.sample, chrom, pos, anc, der, mclass, dp, n_anc, n_der, n_other,
            fwd, rev, n_reads, pct_mq0, n_mq60, call]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
