#!/usr/bin/env python3
"""Place samples on the local tree from read-level pileup tables.

Input is one or more `*_defining.tsv` / `*_control.tsv` outputs of
y_markers_pileup.py -- the per-sample, per-marker tables in results/. Columns are
found by name, so the block files can be passed in any combination and order; a
sample's evidence is whatever rows mention it across all of them.

This is the read-level counterpart to y_path_rank.py --tree. That one scores a
VCF-derived marker_status table and also ranks labels; this one takes the pileup
tables the recent work actually produces and answers only the placement question:
which node is this man at, which nodes is he demonstrably not at, and does the
answer depend on a node no published tree carries.

What it does NOT do is decide anything about power. `ylib.site_call()` already
did that upstream -- one ancestral read is `low_power_1read_ancestral` and is
carried here as a nocall, not as evidence of absence. This script never looks at
read counts, only at calls, so it cannot quietly relax that rule.

Example:
  python3 annotate/y_tree_place.py \\
      --pileup iceman-y/results/testC/testC_*.tsv \\
      --out iceman-y/results/testC/testC_placement.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ytree

SUMMARY_COLUMNS = (
    "sample", "terminal", "terminal_status", "def_derived", "def_derived_transversion",
    "caveats", "excluded", "conflicts", "unpublished_on_path", "path",
)

REQUIRED = ("marker", "call")


def read_pileup(path: Path, scorers: Dict[str, ytree.TreeScorer],
                tree: ytree.Tree, require_qc: bool,
                seen: Dict[tuple, str], disagreements: List[str],
                sample_override: str = "") -> int:
    """Add one pileup file's rows. `seen` deduplicates across files.

    The block files overlap on purpose -- markers/L166_defining.txt and
    markers/yfull_L166_defining.txt share ten positions -- so passing both, which
    is the normal thing to do, would otherwise count the same read twice and
    inflate every derived and ancestral total in the output.
    """
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return 0
        idx = {name: i for i, name in enumerate(h.strip() for h in header)}
        missing = [c for c in REQUIRED if c not in idx]
        if missing:
            # Skip rather than abort: `--pileup iceman-y/results/testC/testC_*.tsv` is the
            # natural invocation and those directories also hold coverage,
            # verdict and params tables. Reported so a renamed column in a real
            # pileup cannot vanish quietly.
            print(f"  skipped {path.name}: not a pileup table (no "
                  f"{', '.join(missing)} column)", file=sys.stderr)
            return 0
        sample_i = idx.get("sample")
        if sample_i is None and not sample_override:
            print(f"  skipped {path.name}: no sample column (pass --sample NAME "
                  "for a single-sample table)", file=sys.stderr)
            return 0
        mclass_i = idx.get("mutation_class")
        qc_i = idx.get("site_qc")
        if require_qc and qc_i is None:
            raise SystemExit(f"ERROR: {path}: --require-site-qc given but no "
                             "site_qc column in this file")

        n = 0
        for row in reader:
            if len(row) <= max(idx[c] for c in REQUIRED):
                continue
            sample = sample_override or row[sample_i].strip()
            marker = row[idx["marker"]].strip()
            if not sample or not marker:
                continue
            call = row[idx["call"]].strip()
            if require_qc and row[qc_i].strip() != "pass":
                call = "no_coverage"

            key = (sample, marker)
            if key in seen:
                if seen[key] != call:
                    disagreements.append(
                        f"{sample} {marker}: {seen[key]} vs {call} ({path.name})")
                continue
            seen[key] = call

            scorer = scorers.get(sample)
            if scorer is None:
                scorer = scorers[sample] = ytree.TreeScorer(tree, name=sample)
            scorer.add_call(marker, call,
                            row[mclass_i] if mclass_i is not None else "")
            n += 1
    return n


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Place samples on the local Y tree from pileup tables.")
    ap.add_argument("--pileup", nargs="+", required=True,
                    help="pileup TSVs from y_markers_pileup.py")
    ap.add_argument("--tree", default="markers/tree_local.tsv", help="tree TSV")
    ap.add_argument("--out", required=True, help="per-node output TSV")
    ap.add_argument("--summary-out",
                    help="one row per sample (default: <out> with .summary.tsv suffix)")
    ap.add_argument("--require-site-qc", action="store_true",
                    help="treat any site not at site_qc=pass as uncovered")
    ap.add_argument("--sample",
                    help="name to use for tables with no sample column, e.g. the "
                         "single-sample Iceman and pooled evidence tables")
    ap.add_argument("--quiet", action="store_true", help="suppress the per-sample report")
    args = ap.parse_args(argv)

    try:
        tree = ytree.load_tree(args.tree)
    except ytree.TreeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    scorers: Dict[str, ytree.TreeScorer] = {}
    seen: Dict[tuple, str] = {}
    disagreements: List[str] = []
    total = 0
    for pattern in args.pileup:
        path = Path(pattern)
        if not path.exists():
            print(f"ERROR: no such file: {path}", file=sys.stderr)
            return 1
        total += read_pileup(path, scorers, tree, args.require_site_qc,
                             seen, disagreements, args.sample or "")

    if not scorers:
        print("ERROR: no sample rows found in the given pileups", file=sys.stderr)
        return 1

    out_path = Path(args.out)
    summary_path = (Path(args.summary_out) if args.summary_out
                    else out_path.with_suffix(out_path.suffix + ".summary.tsv"))

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(ytree.NODE_COLUMNS)
        for name in sorted(scorers):
            w.writerows(scorers[name].rows())

    with summary_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SUMMARY_COLUMNS)
        for name in sorted(scorers):
            s = scorers[name]
            terminal = s.terminal()
            t = terminal[0] if len(terminal) == 1 else ""
            unpublished = [n for n in (s.tree.path_nodes(t) if t else [])
                           if tree.nodes[n].status != "published"]
            w.writerow([
                name,
                ",".join(terminal) if terminal else "-",
                tree.nodes[t].status if t else ("ambiguous" if terminal else "-"),
                s.def_derived[t] if t else 0,
                s.def_derived_tv[t] if t else 0,
                ";".join(s.caveats(t)) if t else "-",
                ",".join(s.excluded()) or "-",
                ",".join(s.conflicts()) or "-",
                ",".join(unpublished) or "-",
                s.ladder(t) if t else "-",
            ])

    if not args.quiet:
        print(f"{total} marker rows, {len(scorers)} samples, tree {tree.path}")
        width = max(len(n) for n in scorers)
        for name in sorted(scorers):
            s = scorers[name]
            line = f"  {name:<{width}}  {s.summary()}"
            if s.excluded():
                line += f"  | not at {', '.join(s.excluded())}"
            if s.conflicts():
                line += f"  | CONFLICT {', '.join(s.conflicts())}"
            print(line)
        if disagreements:
            print(f"  WARNING: {len(disagreements)} marker(s) called differently in "
                  "different input files; the first call was kept:")
            for d in disagreements[:10]:
                print(f"    {d}")
        unmapped = sorted(set().union(*(s.unmapped for s in scorers.values())))
        if unmapped:
            print(f"  {len(unmapped)} marker(s) not in the tree, ignored: "
                  + ", ".join(unmapped[:12]) + (" ..." if len(unmapped) > 12 else ""))
        print(f"  ({ytree.STATUS_LEGEND})")

    print(f"Done. wrote {out_path} and {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
