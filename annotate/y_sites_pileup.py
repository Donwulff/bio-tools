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

from ylib import (detect_mq_ceiling, mapq_audit, mutation_class, region_flag,
                  site_call, site_qc)

INDEL = re.compile(r"[+-](\d+)")


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

    ceiling = detect_mq_ceiling(a.bam)

    cols = ["sample", "chrom", "pos", "anc", "der", "mut_class", "dp", "n_anc",
            "n_der", "n_other", "fwd", "rev", "n_reads", "pct_mq0", "mq_top",
            "n_mq_top", "call", "site_qc", "region"]
    print("\t".join(cols))

    for ln in open(a.sites):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        chrom, pos, anc, der = f[0], int(f[1]), f[2].upper(), f[3].upper()
        # Class is derived from the alleles, not trusted from the input file:
        # the call rules turn on it, so it must not vary with how a site list
        # happened to be written.
        mclass = mutation_class(anc, der)

        _, dp, counts = pileup(a.bam, a.ref, chrom, pos,
                               a.min_mq, a.min_bq, a.max_depth)
        n_reads, n_mq0, n_top = mapq_audit(a.bam, chrom, pos, ceiling)

        n_anc = counts[(anc, "+")] + counts[(anc, "-")]
        n_der = counts[(der, "+")] + counts[(der, "-")]
        n_other = dp - n_anc - n_der
        fwd = sum(v for (b, s), v in counts.items() if s == "+")
        rev = sum(v for (b, s), v in counts.items() if s == "-")

        call = site_call(dp, n_anc, n_der, mclass)
        pct = 100.0 * n_mq0 / n_reads if n_reads else None
        pct_mq0 = f"{pct:.0f}%" if pct is not None else "NA"
        print("\t".join(str(x) for x in [
            a.sample, chrom, pos, anc, der, mclass, dp, n_anc, n_der, n_other,
            fwd, rev, n_reads, pct_mq0, ceiling, n_top, call,
            site_qc(pct, n_top, fwd, rev), region_flag(chrom, pos)]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
