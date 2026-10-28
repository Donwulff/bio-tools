#!/usr/bin/env python3
"""
Compare marker-status TSV outputs from y_haplo_from_markers.py across datasets.

Input rows are expected in this column order:
  CHROM POS ID REF ALT AA HG ISOGG GT GQ DP SOURCE STATUS
Newer tables append REASON; extra columns are accepted. Filtered states remain
distinct from missing records and explicit no-calls.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from marker_schema import CALLABLE, STATUSES, derived_allele


Key = Tuple[str, ...]
Rows = Dict[Key, List[str]]
PAIR_COLUMNS = (
    "dataset_a", "dataset_b", "keys_global_union", "keys_pair_union", "keys_shared",
    "same_status", "different_status", "derived_ancestral_flips", "callable_both",
    "callable_same", "callable_different", "callable_agreement_pct",
    "callable_a_only", "callable_b_only", "callable_unique_locus_pairs",
    "reverse_complemented_keys", "unverified_allele_keys",
)


def parse_input_arg(raw: str) -> Tuple[str, Path]:
    if "=" not in raw:
        raise ValueError(f"--input must be name=path, got: {raw}")
    name, path = raw.split("=", 1)
    name = name.strip()
    if not name or any(c in name for c in "\t\r\n") or not path.strip():
        raise ValueError("--input requires a nonempty dataset name and file path")
    return name, Path(path.strip())


def load_id_filter(path: Path) -> set[str]:
    out: set[str] = set()
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            t = line.strip()
            if not t or t.startswith("#"):
                continue
            out.add(t.split()[0])
    return out


def allele_signature(row: List[str]) -> Optional[Tuple[str, str]]:
    anc = row[5].upper()
    der = derived_allele(row[3].upper(), row[4].upper(), anc)
    if anc == der or not anc or not der or any(c not in "ACGT" for c in anc + der):
        return None
    return anc, der


def load_rows(path: Path, key_mode: str, duplicate_policy: str,
              id_filter: Optional[set[str]], id_re: Optional[re.Pattern]) -> Tuple[Rows, int]:
    rows: Rows = {}
    locus_calls = {}
    collapsed = 0
    with path.open(encoding="utf-8", newline="") as fh:
        for line, row in enumerate(csv.reader(fh, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            if [v.lower() for v in row[:3]] == ["chrom", "pos", "id"]:
                continue
            if len(row) < 13:
                raise ValueError(f"{path}:{line}: expected at least 13 columns")
            marker_id = row[2]
            if id_filter is not None and marker_id not in id_filter:
                continue
            if id_re is not None and id_re.search(marker_id) is None:
                continue
            try:
                pos = int(row[1])
            except ValueError:
                raise ValueError(f"{path}:{line}: invalid marker position") from None
            if pos < 1 or not row[0] or row[12] not in STATUSES:
                raise ValueError(f"{path}:{line}: invalid chromosome, position or status")
            row[1] = str(pos)
            if key_mode == "id" and marker_id in {"", "."}:
                raise ValueError(f"{path}:{line}: ID mode requires a nonempty marker ID")
            key = (marker_id,) if key_mode == "id" else (row[0], row[1], marker_id)
            if key in rows:
                if duplicate_policy == "collapse-identical" and row == rows[key]:
                    collapsed += 1
                    continue
                raise ValueError(f"{path}:{line}: duplicate comparison key {key!r}")
            signature = allele_signature(row)
            if row[12] in CALLABLE:
                if signature is None:
                    raise ValueError(f"{path}:{line}: callable marker has unverifiable ancestral/derived alleles")
                locus = (row[0], row[1], signature)
                if locus in locus_calls and locus_calls[locus] != row[12]:
                    raise ValueError(f"{path}:{line}: conflicting alias calls at the same locus")
                locus_calls[locus] = row[12]
            rows[key] = row
    return rows, collapsed


def compare_pair(a: str, b: str, ra: Rows, rb: Rows, keys: set[Key],
                 allow_complement: bool) -> Dict[str, object]:
    counts: Counter = Counter()
    locus_pairs = set()
    complement = str.maketrans("ACGT", "TGCA")
    for key in keys:
        left, right = ra.get(key), rb.get(key)
        sa, sb = left[12] if left else "missing", right[12] if right else "missing"
        counts["same_status" if sa == sb else "different_status"] += 1
        if {sa, sb} == CALLABLE:
            counts["derived_ancestral_flips"] += 1
        if left is not None or right is not None:
            counts["keys_pair_union"] += 1
        if left is not None and right is not None:
            counts["keys_shared"] += 1
            x, y = allele_signature(left), allele_signature(right)
            if x is None or y is None:
                counts["unverified_allele_keys"] += 1
            elif x != y:
                if allow_complement and all(len(v) == 1 for v in (*x, *y)) and x == tuple(v.translate(complement) for v in y):
                    counts["reverse_complemented_keys"] += 1
                else:
                    raise ValueError(f"incompatible marker alleles at {key!r} between {a!r} and {b!r}")
        ca, cb = sa in CALLABLE, sb in CALLABLE
        if ca and cb:
            counts["callable_both"] += 1
            counts["callable_same" if sa == sb else "callable_different"] += 1
            locus_pairs.add((left[0], left[1], right[0], right[1]))
        elif ca:
            counts["callable_a_only"] += 1
        elif cb:
            counts["callable_b_only"] += 1
    result: Dict[str, object] = {name: counts[name] for name in PAIR_COLUMNS}
    result.update(dataset_a=a, dataset_b=b, keys_global_union=len(keys),
                  callable_unique_locus_pairs=len(locus_pairs),
                  callable_agreement_pct=(f"{100 * counts['callable_same'] / counts['callable_both']:.6f}"
                                          if counts["callable_both"] else "NA"))
    return result


def load_expected(path: Path) -> Dict[str, str]:
    expected: Dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for line, row in enumerate(csv.reader(fh, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            if [x.lower() for x in row] in (["id", "status"], ["id", "expected"]):
                continue
            if len(row) != 2 or row[0] in {"", "."} or row[1] not in CALLABLE:
                raise ValueError(f"{path}:{line}: expected ID and derived/ancestral, separated by a tab")
            if row[0] in expected:
                raise ValueError(f"{path}:{line}: duplicate expected marker ID")
            expected[row[0]] = row[1]
    if not expected:
        raise ValueError(f"{path}: empty expected-state panel")
    return expected


def expected_report(rows_by_ds: Dict[str, Rows], expected: Dict[str, str]):
    details, summaries = [], []
    for name, rows in rows_by_ds.items():
        by_id = {}
        for row in rows.values():
            if row[2] not in expected:
                continue
            if row[2] in by_id:
                raise ValueError(f"{name}: expected marker {row[2]!r} has multiple loci")
            by_id[row[2]] = row
        counts: Counter = Counter()
        for marker_id, wanted in expected.items():
            row = by_id.get(marker_id)
            status = row[12] if row else "missing"
            outcome = ("match" if status == wanted else "conflict") if status in CALLABLE else "unresolved"
            counts[outcome] += 1
            details.append([name, marker_id, wanted, status, outcome,
                            *([row[0], row[1], row[8], row[9], row[10], row[13] if len(row) > 13 else "."]
                              if row else [".", ".", ".", ".", ".", "missing_marker"])])
        summaries.append([name, len(expected), counts["match"] + counts["conflict"],
                          counts["match"], counts["conflict"], counts["unresolved"]])
    return details, summaries


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
        help="locusid=(CHROM,POS,ID) on one assembly; id=unique IDs whose cross-assembly identity you have verified",
    )
    ap.add_argument("--duplicate-policy", choices=["error", "collapse-identical"], default="error")
    ap.add_argument("--allow-reverse-complement", action="store_true",
                    help="allow complemented single-base marker definitions in ID mode; does not verify the coordinate mapping")
    ap.add_argument("--expected-states", help="TSV of ID and expected derived/ancestral state for explicit path QC")
    ap.add_argument("--force", action="store_true", help="replace existing reports (never input files)")
    args = ap.parse_args()
    if args.allow_reverse_complement and args.key_mode != "id":
        ap.error("--allow-reverse-complement requires --key-mode id")
    if not args.out_prefix.strip():
        ap.error("--out-prefix must not be empty")

    datasets: List[Tuple[str, Path]] = []
    try:
        for raw in args.input:
            name, p = parse_input_arg(raw)
            if any(name == previous for previous, _ in datasets):
                raise ValueError(f"duplicate dataset name: {name!r}")
            if not p.is_file():
                raise ValueError(f"input file not found: {p}")
            datasets.append((name, p))
        id_filter = load_id_filter(Path(args.id_list)) if args.id_list else None
        id_re = re.compile(args.id_regex) if args.id_regex else None
        rows_by_ds: Dict[str, Rows] = {}
        keys_union: set[Key] = set()
        for name, p in datasets:
            rows, collapsed = load_rows(p, args.key_mode, args.duplicate_policy, id_filter, id_re)
            if collapsed:
                print(f"WARN: collapsed {collapsed} identical duplicate rows in {name}", file=sys.stderr)
            rows_by_ds[name] = rows
            keys_union.update(rows)
        pair_rows = [compare_pair(a, b, rows_by_ds[a], rows_by_ds[b], keys_union, args.allow_reverse_complement)
                     for (a, _), (b, _) in combinations(datasets, 2)]
        expected = load_expected(Path(args.expected_states)) if args.expected_states else {}
        for marker_id in expected:
            if ((id_filter is not None and marker_id not in id_filter) or
                    (id_re is not None and id_re.search(marker_id) is None)):
                raise ValueError(f"expected marker {marker_id!r} is excluded by the ID filter")
        expected_rows, expected_summaries = expected_report(rows_by_ds, expected) if expected else ([], [])

        outputs = {kind: Path(f"{args.out_prefix}.{kind}.tsv") for kind in
                   (["matrix", "summary", "pairwise", "expected"] if expected else ["matrix", "summary", "pairwise"])}
        if not expected and Path(f"{args.out_prefix}.expected.tsv").exists():
            raise ValueError("an old expected-state report exists; use a fresh prefix or supply --expected-states")
        inputs = {p.resolve() for _, p in datasets}
        inputs.update(Path(p).resolve() for p in (args.id_list, args.expected_states) if p)
        for path in outputs.values():
            if path.resolve() in inputs:
                raise ValueError(f"report would overwrite an input file: {path}")
            if path.exists() and not path.is_file():
                raise ValueError(f"report path is not a file: {path}")
            if (path.exists() or path.is_symlink()) and not args.force:
                raise ValueError(f"output exists: {path}; use a fresh prefix or --force")
        outputs["matrix"].parent.mkdir(parents=True, exist_ok=True)
    except (ValueError, OSError, re.error) as exc:
        raise SystemExit(f"ERROR: {exc}") from None

    matrix_path, summary_path = outputs["matrix"], outputs["summary"]

    # Matrix output
    header = ["CHROM", "POS", "ID"]
    for name, _ in datasets:
        header.extend([f"{name}.status", f"{name}.GT", f"{name}.DP", f"{name}.GQ", f"{name}.HG"])
    for name, _ in datasets:
        header.extend(f"{name}.{field}" for field in ["CHROM", "POS", "REF", "ALT", "AA", "SOURCE", "reason"])

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
            for name, _ in datasets:
                src = rows_by_ds[name].get(key)
                row.extend([src[i] for i in (0, 1, 3, 4, 5, 11)] + [src[13] if len(src) > 13 else "."]
                           if src else ["."] * 7)
            w.writerow(row)

    # Per-dataset status counts
    status_counts: Dict[str, Counter] = {}
    for name, _ in datasets:
        c = Counter()
        for key in keys_union:
            src = rows_by_ds[name].get(key)
            c[src[12] if src else "missing"] += 1
        status_counts[name] = c

    with summary_path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["dataset", "status", "count"])
        for name, _ in datasets:
            for status, n in sorted(status_counts[name].items()):
                w.writerow([name, status, n])
        if expected:
            w.writerow([])
            w.writerow(["dataset", "expected_total", "expected_callable", "expected_match", "expected_conflict", "expected_unresolved"])
            w.writerows(expected_summaries)
        w.writerow([])
        w.writerow(["dataset_a", "dataset_b", "same_status", "different_status", "derived_ancestral_flips"])
        for r in pair_rows:
            w.writerow([r[k] for k in ("dataset_a", "dataset_b", "same_status", "different_status", "derived_ancestral_flips")])

    with outputs["pairwise"].open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=PAIR_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_rows)
    if expected:
        with outputs["expected"].open("w", encoding="utf-8", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["dataset", "ID", "expected", "status", "comparison", "CHROM", "POS", "GT", "GQ", "DP", "reason"])
            writer.writerows(expected_rows)

    print("Done. " + " ".join(f"{kind}={path}" for kind, path in outputs.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
