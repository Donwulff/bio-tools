"""Shared call rules and MAPQ handling for the Y pileup tools.

y_markers_pileup.py and y_sites_pileup.py must make identical calls from
identical evidence -- one resolves markers by name, the other by coordinate,
but a site is a site. Keeping the logic here stops them drifting apart, which
has already happened once: y_markers_pileup.py spent several commits calling
"ancestral" off a single read while y_sites_pileup.py correctly called it
low_power.

Imported by scripts in this directory; Python puts the script's own directory
on sys.path, so a plain `import ylib` works regardless of the caller's cwd.
"""
from __future__ import annotations

import subprocess

DEAMINATION = {("C", "T"), ("G", "A")}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

# Aligner-dependent MAPQ ceiling. This is not cosmetic: bwa aln caps MAPQ at 37
# while bwa mem uses 60, so a "reads at MAPQ 60" test silently evaluates to zero
# for every read of every bwa aln alignment. That defect flagged 100% of sites
# in the Swiss capture data as MARGINAL until it was caught, because the
# threshold had been calibrated on the Iceman's bwa mem BAM.
MQ_CEILING_MEM = 60
MQ_CEILING_ALN = 37


def mutation_class(anc: str, der: str) -> str:
    if len(anc) != 1 or len(der) != 1:
        return "indel/other"
    if (anc, der) in DEAMINATION:
        return "transition(deamination-prone)"
    if (anc, der) in TRANSITIONS:
        return "transition"
    return "transversion"


def detect_mq_ceiling(bam: str, default: int = MQ_CEILING_MEM) -> int:
    """Read the aligner out of @PG and return its maximum MAPQ.

    Detection rather than a flag, because the value is a property of how the
    BAM was made and getting it wrong is silent -- nothing errors, sites just
    quietly fail a quality test they were never able to pass.
    """
    try:
        hdr = subprocess.run(["samtools", "view", "-H", bam],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return default
    for ln in hdr.splitlines():
        if not ln.startswith("@PG"):
            continue
        cl = ln.lower()
        if "bwa mem" in cl:
            return MQ_CEILING_MEM
        if "bwa samse" in cl or "bwa sampe" in cl or "bwa aln" in cl:
            return MQ_CEILING_ALN
    return default


def mapq_audit(bam: str, chrom: str, pos: int, ceiling: int) -> tuple[int, int, int]:
    """(n_reads, n_mq0, n_at_ceiling) for reads overlapping pos, unfiltered.

    Deliberately unfiltered: the point is to see the reads that the -q 25
    pileup threw away. A site with unremarkable depth and 80% MAPQ-0 reads is a
    collapsed repeat, and depth alone cannot reveal that.
    """
    p = subprocess.run(["samtools", "view", bam, f"{chrom}:{pos}-{pos}"],
                       capture_output=True, text=True, check=True)
    n = n0 = ntop = 0
    for ln in p.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 5:
            continue
        if int(f[1]) & 0x900:          # secondary / supplementary
            continue
        mq = int(f[4])
        n += 1
        if mq == 0:
            n0 += 1
        if mq >= ceiling:
            ntop += 1
    return n, n0, ntop


def site_call(dp: int, a: int, d: int, mclass: str) -> str:
    """Pre-registered allele call.

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


# Regions where collapsed repeats produce well-covered but untrustworthy
# pileups. Advisory: the prereg requires these be "flagged, not silently
# included", and a site here can still be real.
NOGO_REGIONS = [
    ("chrY", 11_100_000, 11_700_000, "11.1-11.7Mb"),
    ("chrY", 26_600_000, 26_700_000, "~26.6Mb"),
    ("chrY", 56_690_000, 56_880_000, "56.69-56.88Mb_Yq_het/PAR2"),
]

MAX_PCT_MQ0 = 30.0


def region_flag(chrom: str, pos: int) -> str:
    for c, lo, hi, name in NOGO_REGIONS:
        if chrom == c and lo <= pos <= hi:
            return f"nogo({name})"
    return "ok"


def site_qc(pct_mq0: float | None, n_top: int, fwd: int, rev: int) -> str:
    """Site usability, kept separate from the allele call so that filtering a
    site never destroys the evidence that was filtered.

    The 30% MQ0 threshold was fixed during the Iceman analysis, before any read
    of any later dataset was examined: usable sites topped out at 14% MQ0 and
    rejected ones began at 39%, so any cut in that gap reproduces the same
    split and 30% is the round number in it.

    `n_top` counts reads at the aligner's own MAPQ ceiling, not a hard-coded 60
    -- see detect_mq_ceiling().
    """
    if pct_mq0 is None:
        return "nocall_noreads"
    if pct_mq0 >= MAX_PCT_MQ0:
        return f"REJECT_mapq({pct_mq0:.0f}%_MQ0)"
    if n_top == 0:
        return "MARGINAL_no_unique_reads"
    if fwd == 0 or rev == 0:
        return "MARGINAL_single_strand"
    return "pass"
