#!/usr/bin/env python3
"""Ask which markers of a block co-segregate with a candidate splitting SNP.

The question
------------
A haplotree lists many SNPs as *equivalent* -- defining one node, unorderable
relative to each other. That list is only as good as the samples used to build
it: an equivalence block can only be split by an individual who sits on the
intermediate branch, so a branch with no surviving descendants leaves the block
merged forever in a tree built from modern men.

If an ancient individual does split the block, exactly one thing distinguishes
"a real node was found" from "one SNP mutated twice": whether **other** markers
of the block move with it. A node supported by a single SNP and a recurrence are
the same observation.

What this reports
-----------------
For each marker, the calls in two cohorts that a candidate node predicts should
differ, and a classification:

* ``splits_with_<ref>`` -- derived in the upstream cohort and ancestral in the
  outgroup, the same pattern as the candidate SNP. Every such marker is one more
  SNP the node rests on.
* ``stays_with_block``  -- ancestral in the upstream cohort, i.e. it remains
  below the node and the block is genuinely split rather than relabelled.
* ``uninformative``     -- no coverage, single damage-prone read, or mixed. The
  common outcome at ancient depths and reported as such, never merged into
  either of the above.

An ``uninformative`` marker has not been tested. That distinction is the whole
point: "no other marker co-segregates" is a finding only over markers that were
covered, and this tool separates the two so the count cannot be inflated by
positions nobody could read.

Independence
------------
A documented patriline is one observation, not N. Pool it with
``annotate/y_pool_family.py`` first and pass the pooled table alongside the
per-sample one; the pooled row enters the cohort under its pool name and the
members are excluded by naming the pool rather than the men.

Usage
-----
    annotate/y_cosegregation.py \\
        --calls results/testB/swiss_yfull_L166_defining.tsv \\
        --calls results/testB/swiss_yfull_pooled.tsv \\
        --calls results/testB_unhedged/unhedged_yfull_L166_defining.tsv \\
        --upstream MX210,MX213,SX10,FamilyA \\
        --outgroup I14677,I14678 \\
        --ref-marker Z6219 \\
        --also I5118
"""
from __future__ import annotations

import argparse
import sys

# Calls that assert an allele. Everything else -- no_coverage, the single-read
# nocalls, mixed -- is uninformative by the registered rules and stays that way
# here; this script classifies patterns, it does not re-adjudicate calls.
DERIVED_CALLS = {"DERIVED", "DERIVED_1read_transversion"}
ANCESTRAL_CALLS = {"ancestral"}


def load(paths: list[str]) -> dict[tuple[str, str], dict]:
    """Read per-sample and pooled genotype tables into one {(sample, marker)} map.

    The two formats are told apart by their header rather than by filename:
    y_markers_pileup.py emits `sample/block/marker/...`, y_pool_family.py emits
    `pool/marker/...pooled_call`. Guessing from the path would break silently the
    first time a table is renamed.
    """
    out: dict[tuple[str, str], dict] = {}
    for p in paths:
        with open(p) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {k: i for i, k in enumerate(header)}
            pooled = "pooled_call" in col
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) < len(header):
                    continue
                if pooled:
                    key = (f[col["pool"]], f[col["marker"]])
                    rec = {
                        "call": f[col["pooled_call"]],
                        "anc": f[col["pooled_anc"]],
                        "der": f[col["pooled_der"]],
                        "qc": "pooled:" + f[col["admissible"]],
                        "pos": f[col["pos"]],
                        "class": f[col["mutation_class"]],
                    }
                    # An inadmissible pooled site is a voided sum, not a call.
                    if f[col["admissible"]] != "yes":
                        rec["call"] = "voided_inadmissible"
                else:
                    key = (f[col["sample"]], f[col["marker"]])
                    rec = {
                        "call": f[col["call"]],
                        "anc": f[col["anc_reads"]],
                        "der": f[col["der_reads"]],
                        "qc": f[col["site_qc"]],
                        "pos": f[col["pos"]],
                        "class": f[col["mutation_class"]],
                    }
                out[key] = rec
    return out


def summarise(calls: dict, group: list[str], marker: str) -> tuple[str, str]:
    """(verdict, detail) for one cohort at one marker."""
    der = [s for s in group
           if calls.get((s, marker), {}).get("call") in DERIVED_CALLS]
    anc = [s for s in group
           if calls.get((s, marker), {}).get("call") in ANCESTRAL_CALLS]
    detail = ",".join(
        f"{s}:{calls.get((s, marker), {}).get('anc', '-')}a/"
        f"{calls.get((s, marker), {}).get('der', '-')}d"
        for s in group)
    if der and anc:
        return "split", detail          # cohort disagrees with itself
    if der:
        return "derived", detail
    if anc:
        return "ancestral", detail
    return "nocall", detail


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", action="append", required=True,
                    help="genotype table; repeatable, per-sample or pooled")
    ap.add_argument("--upstream", required=True,
                    help="comma-separated cohort expected DERIVED at the "
                         "candidate node (pool names allowed)")
    ap.add_argument("--outgroup", required=True,
                    help="comma-separated cohort expected ANCESTRAL at it")
    ap.add_argument("--ref-marker", required=True,
                    help="the candidate splitting SNP itself; reported first "
                         "and used as the pattern to match")
    ap.add_argument("--markers", default=None,
                    help="marker name file (default: every marker in --calls)")
    ap.add_argument("--also", default="",
                    help="extra samples reported as context columns only; "
                         "they take no part in the classification")
    ap.add_argument("--out", default="-")
    a = ap.parse_args()

    calls = load(a.calls)
    up = [x for x in a.upstream.split(",") if x]
    out_g = [x for x in a.outgroup.split(",") if x]
    extra = [x for x in a.also.split(",") if x]

    if a.markers:
        names = []
        with open(a.markers) as fh:
            for ln in fh:
                ln = ln.split("#", 1)[0].strip()
                if ln:
                    names.append(ln)
    else:
        seen: set[str] = set()
        names = [m for (_s, m) in calls
                 if not (m in seen or seen.add(m))]

    # The reference marker leads the table: every other row is read against it.
    names = ([a.ref_marker] + [m for m in names if m != a.ref_marker]
             if a.ref_marker in names else names)

    fh_out = sys.stdout if a.out == "-" else open(a.out, "w")
    cols = (["marker", "pos", "mutation_class", "upstream", "outgroup",
             "pattern"]
            + [f"{s}_call" for s in extra]
            + ["upstream_detail", "outgroup_detail"])
    print("\t".join(cols), file=fh_out)

    tally: dict[str, int] = {}
    for m in names:
        any_rec = next((calls[(s, m)] for s in up + out_g + extra
                        if (s, m) in calls), None)
        if any_rec is None:
            pos, mclass = ".", "."
        else:
            pos, mclass = any_rec["pos"], any_rec["class"]

        u, u_det = summarise(calls, up, m)
        o, o_det = summarise(calls, out_g, m)

        if u == "derived" and o == "ancestral":
            pattern = f"splits_with_{a.ref_marker}"
        elif u == "ancestral" and o == "ancestral":
            pattern = "stays_with_block"
        elif u == "split":
            pattern = "cohort_disagrees"
        else:
            pattern = "uninformative"
        tally[pattern] = tally.get(pattern, 0) + 1

        row = [m, pos, mclass, u, o, pattern]
        row += [calls.get((s, m), {}).get("call", "no_row") for s in extra]
        row += [u_det, o_det]
        print("\t".join(row), file=fh_out)

    if fh_out is not sys.stdout:
        fh_out.close()
    print("# " + ", ".join(f"{k}={v}" for k, v in sorted(tally.items())),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
