#!/usr/bin/env python3
"""Apply the registered H1/H2/H0 decision to per-sample genotype tables.

Reads the tables written by annotate/y_genotype_batch.sh and emits one verdict
row per sample, following PREREG_swiss_neolithic_L166.md. The point of doing
this in code is that the H0/H2 distinction survives contact with a results
table:

    H0 is a permitted and expected outcome. It is to be reported as H0, not as
    H2. Absence of derived calls is not an ancestral call. This distinction is
    the single most important commitment in this document.

A human summarising a sparse table is exactly where that collapses -- eight
blank rows and one lone ancestral read reads as "ancestral" unless something
refuses to let it. This script refuses.

Usage:
  annotate/y_prereg_verdict.py --dir results/swiss [--coverage results/swiss/swiss_coverage.tsv]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

# Call strings that constitute positive evidence, from site_call() in
# y_markers_pileup.py / y_sites_pileup.py. Anything not listed here is
# non-informative by construction: no_coverage, low_power_1read_ancestral,
# nocall_damage_prone_1read, nocall_1read_transition.
DERIVED_CALLS = {"DERIVED", "DERIVED_1read_transversion"}
ANCESTRAL_CALLS = {"ancestral"}
MIXED_CALLS = {"mixed", "MIXED_check_paralogy"}


def read_table(path: str) -> list[dict]:
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tally(rows: list[dict], sample: str) -> dict:
    t = defaultdict(int)
    for r in rows:
        if r.get("sample") != sample:
            continue
        c = r.get("call", "")
        t["n"] += 1
        if c in DERIVED_CALLS:
            t["der"] += 1
            t["der_markers"] = t.get("der_markers", "")
        elif c in ANCESTRAL_CALLS:
            t["anc"] += 1
        elif c in MIXED_CALLS:
            t["mixed"] += 1
        elif c == "no_coverage":
            t["nocov"] += 1
        else:
            t["lowpower"] += 1
    return t


def verdict(l166: dict) -> tuple[str, str]:
    """Registered primary decision over the L166-defining set."""
    der, anc, mixed = l166["der"], l166["anc"], l166["mixed"]
    if der and anc:
        return "CONFLICT", f"{der} derived and {anc} ancestral calls -- check paralogy"
    if der:
        return "H1_derived", f"{der} derived call(s), 0 ancestral"
    if anc:
        return "H2_ancestral", f"{anc} ancestral call(s), 0 derived"
    if mixed:
        return "H0_no_power", f"only mixed call(s) ({mixed}); no clean evidence"
    return "H0_no_power", (
        f"no informative call: {l166['nocov']} uncovered, "
        f"{l166['lowpower']} below threshold")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="output dir of y_genotype_batch.sh")
    ap.add_argument("--coverage", default=None)
    ap.add_argument("--prefix", default=None,
                    help="table filename prefix; matches PREFIX in "
                         "y_genotype_batch.sh. Default: auto-detect, which "
                         "keeps the legacy 'swiss' tables readable without "
                         "having to remember what produced them.")
    a = ap.parse_args()

    prefix = a.prefix
    if prefix is None:
        for cand in ("y", "swiss"):
            if os.path.exists(os.path.join(a.dir, f"{cand}_L166_defining.tsv")):
                prefix = cand
                break
        else:
            print(f"no *_L166_defining.tsv in {a.dir}; pass --prefix",
                  file=sys.stderr)
            return 1

    l166 = read_table(os.path.join(a.dir, f"{prefix}_L166_defining.tsv"))
    z6494 = read_table(os.path.join(a.dir, f"{prefix}_Z6494_exclusion.tsv"))
    back = read_table(os.path.join(a.dir, f"{prefix}_backbone_control.tsv"))
    if not l166:
        print(f"no {prefix}_L166_defining.tsv in {a.dir}", file=sys.stderr)
        return 1

    cov = {}
    covpath = a.coverage or os.path.join(a.dir, f"{prefix}_coverage.tsv")
    for r in read_table(covpath):
        cov[r["sample"]] = r.get("chrY_DoC_callable", "NA")

    samples = sorted({r["sample"] for r in l166})
    cols = ["sample", "chrY_DoC_callable", "control_PF3239", "control_derived",
            "control_ancestral", "L166_der", "L166_anc", "L166_nocov",
            "L166_lowpower", "Z6494_der", "Z6494_anc", "verdict", "basis"]
    print("\t".join(cols))

    for s in samples:
        tl, tz, tb = tally(l166, s), tally(z6494, s), tally(back, s)
        pf = next((r["call"] for r in back
                   if r.get("sample") == s and r.get("marker") == "PF3239"), "absent")
        v, basis = verdict(tl)

        # Registered: "A sample that fails this control is excluded from
        # interpretation." Failing means contradicting the published backbone --
        # an ancestral call where the publication reports derived. An uncovered
        # control is not a failure, it is no information, and is reported as
        # such rather than used to discard the sample.
        if tb["anc"]:
            v = "EXCLUDED_control_fail"
            basis = (f"ancestral at {tb['anc']} backbone marker(s) where the "
                     f"publication reports derived; not interpretable")
        elif tb["der"] == 0:
            basis += "; backbone control uncovered, chain unverified for this sample"

        # The exclusion set is reported but never converts an H0 into a verdict:
        # being off the Iceman's branch in the other direction is a separate
        # claim from being ancestral at L166.
        if tz["der"]:
            basis += f"; DERIVED at {tz['der']} Z6494 site(s) -- off-branch"

        print("\t".join(str(x) for x in [
            s, cov.get(s, "NA"), pf, tb["der"], tb["anc"],
            tl["der"], tl["anc"], tl["nocov"], tl["lowpower"],
            tz["der"], tz["anc"], v, basis]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
