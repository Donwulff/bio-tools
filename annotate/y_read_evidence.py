#!/usr/bin/env python3
"""Dump the individual reads behind a call, plus the library's damage profile.

Why this exists
---------------
A pileup collapses reads into counts, and counts cannot distinguish "4 molecules
independently carry the derived allele" from "4 molecules are damaged at the same
position". For a C>T or G>A site those are the same tally and opposite
conclusions. The registered rules already refuse a *single* damage-prone read,
but they have nothing to say about four, and at some point a claim has to rest on
looking at the reads.

Two things are needed to interpret a site and this prints both:

* **per-read evidence** -- strand, length, the site's distance from each end of
  the original molecule, base quality, MAPQ, and how many other mismatches the
  same read carries. Deamination is a *terminal* process, so a derived allele
  sitting 20 bp from both ends of every read is not damage, while one sitting at
  position 1 of every read almost certainly is.
* **the library's own damage rate** -- the positional C>T and G>A frequencies
  measured from the same BAM. Without it "these reads are terminal" is
  unquantified: terminal position only matters in proportion to how often this
  library actually deaminates there. A UDG-treated library can put a read at
  position 1 with almost no risk.

Orientation is handled explicitly. A read stored reverse-complemented in the BAM
has its original 5' end at the *right* of the alignment, so naive "distance from
alignment start" mislabels exactly the reads that matter. ``dist_5p``/``dist_3p``
below are distances along the original molecule.

Usage
-----
    annotate/y_read_evidence.py --bam MX210.rmdup.bam --ref hg38.fa \\
        --site chrY:13782251 --anc C --der T

    annotate/y_read_evidence.py --bam MX210.rmdup.bam --ref hg38.fa \\
        --damage-profile --region chrY --max-reads 200000
"""
from __future__ import annotations

import argparse
import collections
import subprocess
import sys

CIGAR_CONSUMES_QUERY = set("MIS=X")
CIGAR_CONSUMES_REF = set("MDN=X")


def parse_cigar(cig: str) -> list[tuple[int, str]]:
    ops, n = [], ""
    for c in cig:
        if c.isdigit():
            n += c
        else:
            ops.append((int(n), c))
            n = ""
    return ops


def query_index_at_ref(pos_ref_start: int, cigar: str, target: int) -> int | None:
    """0-based offset into SEQ of the base aligned to `target` (1-based ref)."""
    q, r = 0, pos_ref_start
    for ln, op in parse_cigar(cigar):
        if op in ("M", "=", "X"):
            if r <= target < r + ln:
                return q + (target - r)
            q += ln
            r += ln
        elif op in ("I", "S"):
            q += ln
        elif op in ("D", "N"):
            if r <= target < r + ln:
                return None          # deleted in this read
            r += ln
        elif op == "H":
            continue
    return None


def samtools_view(bam: str, region: str, extra: list[str] | None = None):
    cmd = ["samtools", "view"] + (extra or []) + [bam, region]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                         stderr=subprocess.DEVNULL)
    assert p.stdout is not None
    for ln in p.stdout:
        yield ln.rstrip("\n").split("\t")
    p.wait()


def ref_base(ref: str, chrom: str, pos: int) -> str:
    p = subprocess.run(["samtools", "faidx", ref, f"{chrom}:{pos}-{pos}"],
                       capture_output=True, text=True)
    return "".join(p.stdout.split("\n")[1:]).upper() or "?"


def per_read(a) -> int:
    chrom, pos = a.site.split(":")
    pos = int(pos)
    rb = ref_base(a.ref, chrom, pos)
    anc = a.anc or rb
    print(f"# {a.site}  ref={rb}  anc={anc}  der={a.der or '?'}  bam={a.bam}",
          file=sys.stderr)

    cols = ["read", "strand", "read_len", "base", "bq", "mapq",
            "dist_5p", "dist_3p", "nm", "other_ct_ga", "verdict"]
    print("\t".join(cols))

    n = collections.Counter()
    for f in samtools_view(a.bam, f"{chrom}:{pos}-{pos}"):
        if len(f) < 11:
            continue
        qname, flag, spos, mapq, cig, seq, qual = (
            f[0], int(f[1]), int(f[3]), int(f[4]), f[5], f[9], f[10])
        if flag & 0x904:                       # unmapped/secondary/supplementary
            continue
        if mapq < a.min_mq:
            continue
        qi = query_index_at_ref(spos, cig, pos)
        if qi is None or qi >= len(seq):
            continue
        base = seq[qi].upper()
        bq = ord(qual[qi]) - 33
        rev = bool(flag & 16)
        L = len(seq)
        # Distances along the ORIGINAL molecule. A reverse-strand read is stored
        # reverse-complemented, so its original 5' end is at the right-hand side.
        d_from_left, d_from_right = qi, L - 1 - qi
        dist_5p = d_from_right if rev else d_from_left
        dist_3p = d_from_left if rev else d_from_right

        nm = next((int(x.split(":")[-1]) for x in f[11:]
                   if x.startswith("NM:i:")), -1)

        # Other deamination-pattern mismatches on the same molecule. A genuinely
        # damaged read usually carries more than one; a read whose only mismatch
        # is the site of interest is weaker evidence for damage.
        other = 0
        md = next((x.split(":", 2)[2] for x in f[11:] if x.startswith("MD:Z:")),
                  None)
        if md is not None:
            other = sum(1 for c in md if c in "CG") - (1 if base != rb else 0)
            other = max(other, 0)

        if base == anc:
            verdict = "ancestral"
        elif a.der and base == a.der:
            dmg = (rb, base) in (("C", "T"), ("G", "A"))
            verdict = f"derived{'(damage-pattern)' if dmg else ''}"
        else:
            verdict = "other"
        n[verdict.split("(")[0]] += 1

        print("\t".join(str(x) for x in [
            qname[:28], "-" if rev else "+", L, base, bq, mapq,
            dist_5p, dist_3p, nm, other, verdict]))

    print(f"# totals: {dict(n)}", file=sys.stderr)
    return 0


def damage_profile(a) -> int:
    """Positional C>T (5') and G>A (3') rates, measured from this BAM's MD tags."""
    ct = collections.Counter()
    ga = collections.Counter()
    cov5 = collections.Counter()
    cov3 = collections.Counter()
    seen = 0
    for f in samtools_view(a.bam, a.region):
        if len(f) < 11 or seen >= a.max_reads:
            break
        flag = int(f[1])
        if flag & 0x904 or int(f[4]) < a.min_mq:
            continue
        seq, rev = f[9], bool(flag & 16)
        L = len(seq)
        md = next((x.split(":", 2)[2] for x in f[11:] if x.startswith("MD:Z:")),
                  None)
        if md is None:
            continue
        seen += 1
        for i in range(min(a.ends, L)):
            cov5[i] += 1
            cov3[i] += 1
        # Reconstruct reference bases at mismatches from MD, walking SEQ.
        qi, num, k = 0, "", 0
        while k < len(md):
            c = md[k]
            if c.isdigit():
                num += c
                k += 1
                continue
            if num:
                qi += int(num)
                num = ""
            if c == "^":
                k += 1
                while k < len(md) and md[k].isalpha():
                    k += 1
                continue
            refb = c.upper()
            if qi < L:
                obs = seq[qi].upper()
                d5 = (L - 1 - qi) if rev else qi
                d3 = qi if rev else (L - 1 - qi)
                if refb == "C" and obs == "T" and d5 < a.ends:
                    ct[d5] += 1
                if refb == "G" and obs == "A" and d3 < a.ends:
                    ga[d3] += 1
            qi += 1
            k += 1

    print(f"# {a.bam} {a.region}: {seen} reads, MAPQ>={a.min_mq}",
          file=sys.stderr)
    print("pos_from_end\tC_to_T_5p\trate_5p\tG_to_A_3p\trate_3p")
    for i in range(a.ends):
        r5 = ct[i] / cov5[i] if cov5[i] else 0.0
        r3 = ga[i] / cov3[i] if cov3[i] else 0.0
        print(f"{i}\t{ct[i]}\t{r5:.4f}\t{ga[i]}\t{r3:.4f}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", default=None)
    ap.add_argument("--site", default=None, help="chrom:pos")
    ap.add_argument("--anc", default=None)
    ap.add_argument("--der", default=None)
    ap.add_argument("--min-mq", type=int, default=25)
    ap.add_argument("--damage-profile", action="store_true")
    ap.add_argument("--region", default="chrY")
    ap.add_argument("--ends", type=int, default=10,
                    help="positions from each end to profile (default 10)")
    ap.add_argument("--max-reads", type=int, default=200000)
    a = ap.parse_args()

    if a.damage_profile:
        return damage_profile(a)
    if not a.site or not a.ref:
        print("--site and --ref required unless --damage-profile", file=sys.stderr)
        return 1
    return per_read(a)


if __name__ == "__main__":
    sys.exit(main())
