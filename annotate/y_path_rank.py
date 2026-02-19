#!/usr/bin/env python3
"""
Rank Y-haplogroup candidates from marker status tables.

Input: *.marker_status.tsv from y_haplo_from_markers.py
Columns:
  1 CHROM
  2 POS
  3 ID
  4 REF
  5 ALT
  6 AA
  7 HG
  8 ISOGG
  9 GT
 10 GQ
 11 DP
 12 SOURCE
 13 STATUS
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Set, Tuple


NOISE_PATTERNS = (
    "unknown",
    "not_listed",
    "approx.",
    "aprrox.",
    "approx_",
    "approx",
)


def is_deamination_transition(ref: str, alt: str) -> bool:
    return (ref == "C" and alt == "T") or (ref == "G" and alt == "A")


def is_transition(ref: str, alt: str) -> bool:
    return (ref, alt) in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def normalize_token(token: str) -> str:
    t = token.strip().strip('"').strip("'")
    t = t.replace("Approx._hg:_", "")
    t = t.replace("Aprrox._hg:_", "")
    t = t.replace("Approx._Hg:_", "")
    t = t.replace("Approx._", "")
    t = t.replace("_", "")
    t = t.replace(" ", "")
    return t


def tokenize_labels(hg: str, isogg: str, clade_prefix: str) -> List[str]:
    raw = ",".join([hg or "", isogg or ""])
    parts = re.split(r"[;,>]", raw)
    out: List[str] = []
    for p in parts:
        t = normalize_token(p)
        if not t or t == ".":
            continue
        low = t.lower()
        if any(np in low for np in NOISE_PATTERNS):
            continue
        if not t.startswith(clade_prefix):
            continue
        out.append(t)
    return out


def hierarchical_prefixes(token: str, clade_prefix: str) -> List[str]:
    """
    Return path-like prefixes for tokens that look hierarchical (e.g. G2a2a1a2a1a1b).
    For label-like SNP names (e.g. G-Z6208), return [token].
    """
    t = token.rstrip("~*")
    prefixes: List[str] = []

    # Case: G2a-L91 -> include base hierarchy from G2a and the full label token.
    m_hy = re.match(rf"^({re.escape(clade_prefix)}[0-9A-Za-z]+)-", t)
    if m_hy:
        base = m_hy.group(1)
        prefixes.extend(hierarchical_prefixes(base, clade_prefix))
        prefixes.append(t)
        return list(dict.fromkeys(prefixes))

    # Pure hierarchical ISOGG-like token.
    m = re.match(rf"^{re.escape(clade_prefix)}([0-9A-Za-z]+)$", t)
    if not m:
        return [t]
    suffix = m.group(1)
    cur = clade_prefix
    prefixes.append(cur)
    for part in re.findall(r"[0-9]+|[A-Za-z]+", suffix):
        cur += part
        prefixes.append(cur)
    return list(dict.fromkeys(prefixes))


def main() -> int:
    ap = argparse.ArgumentParser(description="Rank clade paths from marker status rows.")
    ap.add_argument("--input", required=True, help="marker_status.tsv from y_haplo_from_markers.py")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--clade-prefix", default="G", help="target clade prefix (default: G)")
    ap.add_argument(
        "--deam-derived-weight",
        type=float,
        default=0.5,
        help="weight for derived C>T / G>A markers (default: 0.5)",
    )
    ap.add_argument(
        "--normal-derived-weight",
        type=float,
        default=1.0,
        help="weight for other derived markers (default: 1.0)",
    )
    ap.add_argument(
        "--ancestral-weight",
        type=float,
        default=-1.0,
        help="weight for ancestral markers at candidate labels (default: -1.0)",
    )
    ap.add_argument(
        "--nocall-weight",
        type=float,
        default=0.0,
        help="weight for nocall markers at candidate labels (default: 0.0)",
    )
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"ERROR: input file not found: {in_path}")

    derived_total = defaultdict(int)
    derived_transversion = defaultdict(int)
    derived_deamination = defaultdict(int)
    ancestral_total = defaultdict(int)
    nocall_total = defaultdict(int)
    ambiguous_total = defaultdict(int)
    other_total = defaultdict(int)
    score_total = defaultdict(float)
    candidates_seen: Set[str] = set()

    with in_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 13:
                continue
            ref = row[3]
            alt = row[4].split(",")[0]
            hg = row[6]
            isogg = row[7]
            status = row[12]

            labels = tokenize_labels(hg, isogg, args.clade_prefix)
            if not labels:
                continue

            cand_for_row: Set[str] = set()
            for label in labels:
                for c in hierarchical_prefixes(label, args.clade_prefix):
                    cand_for_row.add(c)
            if not cand_for_row:
                continue

            deam = is_deamination_transition(ref, alt)
            transv = not is_transition(ref, alt)

            for c in cand_for_row:
                candidates_seen.add(c)
                if status == "derived":
                    w = args.deam_derived_weight if deam else args.normal_derived_weight
                    derived_total[c] += 1
                    score_total[c] += w
                    if deam:
                        derived_deamination[c] += 1
                    if transv:
                        derived_transversion[c] += 1
                elif status == "ancestral":
                    ancestral_total[c] += 1
                    score_total[c] += args.ancestral_weight
                elif status == "nocall":
                    nocall_total[c] += 1
                    score_total[c] += args.nocall_weight
                elif status == "ambiguous":
                    ambiguous_total[c] += 1
                else:
                    other_total[c] += 1

    ranked = sorted(
        candidates_seen,
        key=lambda c: (score_total[c], derived_total[c], -ancestral_total[c], c),
        reverse=True,
    )

    out_path = Path(args.out)
    with out_path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(
            [
                "candidate",
                "score",
                "derived_total",
                "derived_transversion",
                "derived_deamination",
                "ancestral_total",
                "nocall_total",
                "ambiguous_total",
                "other_total",
            ]
        )
        for c in ranked:
            w.writerow(
                [
                    c,
                    f"{score_total[c]:.3f}",
                    derived_total[c],
                    derived_transversion[c],
                    derived_deamination[c],
                    ancestral_total[c],
                    nocall_total[c],
                    ambiguous_total[c],
                    other_total[c],
                ]
            )

    print(f"Done. wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
