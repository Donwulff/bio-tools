#!/usr/bin/env python3
"""
Score Y-clade consistency from derived marker tables produced by y_haplo_from_vcf.sh.

Input table format (TSV), expected columns:
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
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


NOISE_PATTERNS = (
    "unknown",
    "not_listed",
    "approx.",
    "aprrox.",
    "approx_",
    "approx",
)


def normalize_token(token: str) -> str:
    t = token.strip().strip('"').strip("'")
    t = t.replace("Approx._hg:_", "")
    t = t.replace("Aprrox._hg:_", "")
    t = t.replace("Approx._Hg:_", "")
    t = t.replace("Approx._hg__", "")
    t = t.replace("Approx._", "")
    t = t.replace("_", "")
    t = t.replace(" ", "")
    t = t.strip()
    return t


def split_tokens(hg: str, isogg: str) -> list[str]:
    raw = ",".join([hg or "", isogg or ""])
    parts = re.split(r"[;,]", raw)
    out: list[str] = []
    for part in parts:
        t = normalize_token(part)
        if not t or t == ".":
            continue
        low = t.lower()
        if any(p in low for p in NOISE_PATTERNS):
            continue
        out.append(t)
    return out


def major_clade(token: str) -> str | None:
    m = re.match(r"^([A-Z])", token)
    return m.group(1) if m else None


def strip_star(clade: str) -> str:
    return clade.replace("*", "")


def starts_with_clade(token: str, clade: str) -> bool:
    # Accept exact or descendant forms, e.g. G2a2 startswith G2a; G-Z6208 startswith G-Z6208.
    return token.startswith(clade)


def token_is_upstream_of_candidate(token: str, clade: str) -> bool:
    # Upstream if token is a prefix of candidate and shares major clade, e.g. G or G2a for G2a2a1...
    if not clade.startswith(token):
        return False
    m1 = major_clade(token)
    m2 = major_clade(clade)
    return m1 is not None and m1 == m2


@dataclass
class Counts:
    total_rows: int = 0
    informative_rows: int = 0
    support_specific: int = 0
    support_upstream: int = 0
    support_any: int = 0
    conflict_major: int = 0
    conflict_exclusive: int = 0
    mixed_support_conflict: int = 0


def score_dataset(path: Path, candidate: str) -> Counts:
    c = Counts()
    cand = strip_star(candidate)
    cand_major = major_clade(cand)

    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 8:
                continue
            c.total_rows += 1
            tokens = split_tokens(row[6], row[7])
            if not tokens:
                continue
            c.informative_rows += 1

            support_specific = any(starts_with_clade(strip_star(t), cand) for t in tokens)
            support_upstream = (not support_specific) and any(
                token_is_upstream_of_candidate(strip_star(t), cand) for t in tokens
            )
            support_any = support_specific or support_upstream

            conflict_major = False
            if cand_major is not None:
                for t in tokens:
                    tm = major_clade(t)
                    if tm is not None and tm != cand_major:
                        conflict_major = True
                        break

            if support_specific:
                c.support_specific += 1
            if support_upstream:
                c.support_upstream += 1
            if support_any:
                c.support_any += 1
            if conflict_major:
                c.conflict_major += 1
                if support_any:
                    c.mixed_support_conflict += 1
                else:
                    c.conflict_exclusive += 1

    return c


def parse_input_arg(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("input must be label=path")
    label, path = value.split("=", 1)
    if not label:
        raise argparse.ArgumentTypeError("input label must not be empty")
    p = Path(path)
    if not p.exists():
        raise argparse.ArgumentTypeError(f"input file not found: {path}")
    return label, p


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Score candidate Y clades against derived marker TSV tables."
    )
    ap.add_argument(
        "--input",
        action="append",
        required=True,
        type=parse_input_arg,
        metavar="LABEL=PATH",
        help="Derived TSV input (repeatable)",
    )
    ap.add_argument(
        "--candidate",
        action="append",
        required=True,
        metavar="CLADE",
        help="Candidate clade (repeatable), e.g. G2a2a1a2a1a1b",
    )
    ap.add_argument(
        "--out",
        default="-",
        metavar="PATH",
        help="Output TSV path (default: stdout)",
    )
    args = ap.parse_args()

    out_fh = sys.stdout if args.out == "-" else open(args.out, "w", encoding="utf-8", newline="")
    try:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(
            [
                "dataset",
                "candidate",
                "total_rows",
                "informative_rows",
                "support_specific",
                "support_upstream",
                "support_any",
                "conflict_major",
                "conflict_exclusive",
                "mixed_support_conflict",
                "net_support_any_minus_conflict_exclusive",
            ]
        )
        for label, path in args.input:
            for cand in args.candidate:
                c = score_dataset(path, cand)
                writer.writerow(
                    [
                        label,
                        cand,
                        c.total_rows,
                        c.informative_rows,
                        c.support_specific,
                        c.support_upstream,
                        c.support_any,
                        c.conflict_major,
                        c.conflict_exclusive,
                        c.mixed_support_conflict,
                        c.support_any - c.conflict_exclusive,
                    ]
                )
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
