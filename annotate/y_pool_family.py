#!/usr/bin/env python3
"""Sum read counts across a documented patrilineal kin group, per marker.

Why this is not the pooling the pre-registration refuses
-------------------------------------------------------
``PREREG_swiss_neolithic_L166.md`` declines to pool reads across the Oberbipp
males, and the reason it gives is specific: combining unrelated men to ask "are
they L166" assumes they share a branch, which is the hypothesis under test.

That reason does not apply to a set *independently documented* to share a Y
chromosome. Male-male first-degree pairs are father-son or full brothers, so
their Y chromosomes are one chromosome sampled repeatedly. Summing their reads
is then the correct estimator rather than a power rescue -- but only for the
question "what is this patriline's terminal haplogroup", and the answer counts
as exactly **one** observation. It can never become a proportion over
individuals. See Amendment 2, section D.

Membership must come from published kinship, never from the genotypes being
pooled; ``markers/family_a_members.txt`` records the source text.

Admissibility
-------------
A pooled site is voided if **any** member with reads at that position fails
``site_qc`` with a REJECT. One member's collapsed repeat contaminates the sum,
and the sum is what gets called, so the strictest contributor governs. Members
with no reads are simply absent from the sum.

Per-sample contributions are printed alongside every pooled row so a reader can
undo the pooling and check it.

Example
-------
    annotate/y_pool_family.py \\
        --genotypes results/swiss15/swiss_L166_defining.tsv \\
        --members markers/family_a_members.txt \\
        --pool-name FamilyA
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict

from ylib import site_call


def load_members(path: str) -> list[str]:
    out: list[str] = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.split("#", 1)[0].strip()
            if ln:
                out.append(ln)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True,
                    help="per-sample TSV from annotate/y_genotype_batch.sh")
    ap.add_argument("--members", required=True,
                    help="file of sample IDs in the kin group, one per line")
    ap.add_argument("--pool-name", default="pool")
    ap.add_argument("--out", default="-")
    args = ap.parse_args()

    members = load_members(args.members)
    mset = set(members)

    with open(args.genotypes) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"no rows in {args.genotypes}")
    for need in ("sample", "marker", "anc_reads", "der_reads", "site_qc"):
        if need not in rows[0]:
            sys.exit(f"{args.genotypes} has no '{need}' column; "
                     f"expected output of annotate/y_genotype_batch.sh")

    seen_samples = {r["sample"] for r in rows}
    absent = [m for m in members if m not in seen_samples]
    if absent:
        print(f"# members with no genotype rows ({len(absent)}): {','.join(absent)}",
              file=sys.stderr)

    by_marker: dict[str, list[dict]] = defaultdict(list)
    order: list[str] = []
    for r in rows:
        if r["sample"] not in mset:
            continue
        if r["marker"] not in by_marker:
            order.append(r["marker"])
        by_marker[r["marker"]].append(r)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    cols = ["pool", "marker", "chrom", "pos", "anc>der", "mutation_class",
            "n_contributors", "pooled_dp", "pooled_anc", "pooled_der",
            "pooled_other", "pooled_call", "admissible", "contributions",
            "region"]
    print("\t".join(cols), file=out)

    for marker in order:
        rs = by_marker[marker]
        ref = rs[0]
        anc = der = other = 0
        contribs: list[str] = []
        vetoes: list[str] = []
        for r in rs:
            a = int(r["anc_reads"] or 0)
            d = int(r["der_reads"] or 0)
            o = int(r["other_reads"] or 0)
            nreads = int(r["n_reads"] or 0) if (r.get("n_reads") or "").isdigit() else 0
            qc = r["site_qc"]
            # A member with reads at the position governs the sum whether or not
            # those reads survived the pileup's MAPQ filter -- reads rejected for
            # mapping quality are evidence the locus is bad, not evidence of
            # nothing.
            if nreads > 0 and qc.startswith("REJECT"):
                vetoes.append(f"{r['sample']}:{qc}")
                continue
            if a or d or o:
                contribs.append(f"{r['sample']}:{a}a/{d}d")
            anc += a
            der += d
            other += o

        dp = anc + der + other
        if vetoes:
            call = "VOID_contributor_rejected"
            admissible = "no(" + ";".join(vetoes) + ")"
        else:
            call = site_call(dp, anc, der, ref["mutation_class"])
            admissible = "yes"

        print("\t".join([
            args.pool_name, marker, ref["chrom"], ref["pos"],
            ref["anc>der"], ref["mutation_class"],
            str(len(contribs)), str(dp), str(anc), str(der), str(other),
            call, admissible,
            ",".join(contribs) if contribs else ".",
            ref.get("region", "."),
        ]), file=out)

    if out is not sys.stdout:
        out.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
