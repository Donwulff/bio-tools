#!/usr/bin/env python3
"""
Compare marker-status TSV outputs from y_haplo_from_markers.py across datasets.

Input rows are expected in this column order:
  CHROM POS ID REF ALT AA HG ISOGG GT GQ DP SOURCE STATUS
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def parse_input_arg(raw: str) -> Tuple[str, Path]:
    if "=" not in raw:
        raise ValueError(f"--input must be name=path, got: {raw}")
    name, path = raw.split("=", 1)
    name = name.strip()
    p = Path(path.strip())
    if not name:
        raise ValueError(f"empty dataset name in --input: {raw}")
    return name, p


def load_id_filter(path: Path) -> set[str]:
    out: set[str] = set()
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            t = line.strip()
            if not t or t.startswith("#"):
                continue
            out.add(t.split()[0])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Compare marker-status tables across datasets.")
    ap.add_argument(
        "--input",
        action="append",
        required=True,
        help="dataset_name=marker_status.tsv (repeat)",
    )
    ap.add_argument("--out-prefix", required=True, help="output file prefix")
    ap.add_argument("--id-list", help="optional file with marker IDs (one per line)")
    ap.add_argument("--id-regex", help="optional regex filter on marker ID")
    ap.add_argument(
        "--key-mode",
        choices=["locusid", "id"],
        default="locusid",
        help="comparison key: locusid=(CHROM,POS,ID) default, or id=(ID only, useful across assemblies)",
    )
    args = ap.parse_args()

    datasets: List[Tuple[str, Path]] = []
    for raw in args.input:
        name, p = parse_input_arg(raw)
        if not p.exists():
            raise SystemExit(f"ERROR: input not found: {p}")
        datasets.append((name, p))

    id_filter = None
    if args.id_list:
        id_filter = load_id_filter(Path(args.id_list))
    id_re = re.compile(args.id_regex) if args.id_regex else None

    # key = (chrom,pos,id) or (id,)
    rows_by_ds: Dict[str, Dict[Tuple[str, str, str], List[str]]] = {}
    keys_union = set()
    for name, p in datasets:
        d: Dict[Tuple[str, str, str], List[str]] = {}
        with p.open("r", encoding="utf-8", newline="") as fh:
            r = csv.reader(fh, delimiter="\t")
            for row in r:
                if len(row) < 13:
                    continue
                marker_id = row[2]
                if id_filter is not None and marker_id not in id_filter:
                    continue
                if id_re is not None and id_re.search(marker_id) is None:
                    continue
                if args.key_mode == "id":
                    key = (marker_id,)
                else:
                    key = (row[0], row[1], marker_id)
                d[key] = row
                keys_union.add(key)
        rows_by_ds[name] = d

    out_prefix = Path(args.out_prefix)
    matrix_path = out_prefix.with_suffix(".matrix.tsv")
    summary_path = out_prefix.with_suffix(".summary.tsv")

    # Matrix output
    header = ["CHROM", "POS", "ID"]
    for name, _ in datasets:
        header.extend([f"{name}.status", f"{name}.GT", f"{name}.DP", f"{name}.GQ", f"{name}.HG"])

    with matrix_path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        if args.key_mode == "id":
            key_iter = sorted(keys_union, key=lambda k: k[0])
        else:
            key_iter = sorted(keys_union, key=lambda k: (k[0], int(k[1]), k[2]))
        for key in key_iter:
            if args.key_mode == "id":
                row = [".", ".", key[0]]
            else:
                row = [key[0], key[1], key[2]]
            for name, _ in datasets:
                src = rows_by_ds[name].get(key)
                if src is None:
                    row.extend(["missing", ".", ".", ".", "."])
                else:
                    row.extend([src[12], src[8], src[10], src[9], src[6]])
            w.writerow(row)

    # Per-dataset status counts
    status_counts: Dict[str, Counter] = {}
    for name, _ in datasets:
        c = Counter()
        for key in keys_union:
            src = rows_by_ds[name].get(key)
            c[src[12] if src else "missing"] += 1
        status_counts[name] = c

    # Pairwise agreement
    pair_rows: List[List[str]] = []
    for i in range(len(datasets)):
        for j in range(i + 1, len(datasets)):
            a = datasets[i][0]
            b = datasets[j][0]
            same = 0
            diff = 0
            flip = 0  # derived<->ancestral flips
            for key in keys_union:
                ra = rows_by_ds[a].get(key)
                rb = rows_by_ds[b].get(key)
                sa = ra[12] if ra else "missing"
                sb = rb[12] if rb else "missing"
                if sa == sb:
                    same += 1
                else:
                    diff += 1
                if {sa, sb} == {"derived", "ancestral"}:
                    flip += 1
            pair_rows.append([a, b, str(same), str(diff), str(flip)])

    with summary_path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["dataset", "status", "count"])
        for name, _ in datasets:
            for status, n in sorted(status_counts[name].items()):
                w.writerow([name, status, n])
        w.writerow([])
        w.writerow(["dataset_a", "dataset_b", "same_status", "different_status", "derived_ancestral_flips"])
        for r in pair_rows:
            w.writerow(r)

    print(f"Done. matrix={matrix_path} summary={summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
