#!/usr/bin/env python3
"""
Summarize insert-size distribution from samtools stats output (IS lines).

Usage examples:
  samtools stats sample.bam | ./util/is_stats.py
  ./util/is_stats.py sample.stats
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize insert-size distribution from samtools stats (IS lines)."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="samtools stats file (default: stdin).",
    )
    parser.add_argument(
        "--pcts",
        default="50,90,99",
        help="Comma-separated percentiles to report (default: 50,90,99).",
    )
    parser.add_argument(
        "--tail",
        type=int,
        default=200,
        help="Tail threshold (bp) to report fraction >= threshold (default: 200).",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=5,
        help="Smoothing half-window for peak detection (default: 5).",
    )
    parser.add_argument(
        "--peak-threshold",
        type=float,
        default=0.10,
        help="Peak height threshold as fraction of max smoothed count (default: 0.10).",
    )
    return parser.parse_args()


def load_is_counts(stream) -> Dict[int, int]:
    counts: Dict[int, int] = {}
    for line in stream:
        if not line.startswith("IS"):
            continue
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        try:
            size = int(parts[1])
            count = int(parts[2])
        except ValueError:
            continue
        counts[size] = counts.get(size, 0) + count
    return counts


def percentiles(counts: Dict[int, int], pcts: List[float]) -> Dict[float, int]:
    total = sum(counts.values())
    if total == 0:
        return {p: 0 for p in pcts}
    targets = {p: total * (p / 100.0) for p in pcts}
    result: Dict[float, int] = {}
    cum = 0
    for size in sorted(counts):
        cum += counts[size]
        for p, t in targets.items():
            if p not in result and cum >= t:
                result[p] = size
    for p in pcts:
        result.setdefault(p, 0)
    return result


def smooth_counts(counts: Dict[int, int], window: int) -> Tuple[List[int], List[int]]:
    if not counts:
        return [], []
    max_size = max(counts)
    arr = [0] * (max_size + 1)
    for size, count in counts.items():
        if size >= 0:
            arr[size] = count
    smoothed = []
    for i in range(len(arr)):
        lo = max(0, i - window)
        hi = min(len(arr), i + window + 1)
        smoothed.append(sum(arr[lo:hi]))
    return arr, smoothed


def detect_peaks(smoothed: List[int], threshold: float) -> List[int]:
    if not smoothed:
        return []
    max_val = max(smoothed)
    if max_val == 0:
        return []
    peaks = []
    for i in range(1, len(smoothed) - 1):
        if smoothed[i - 1] < smoothed[i] > smoothed[i + 1] and smoothed[i] >= max_val * threshold:
            peaks.append(i)
    return peaks


def main() -> int:
    args = parse_args()
    if args.input == "-":
        counts = load_is_counts(sys.stdin)
    else:
        with open(args.input, "r", encoding="utf-8") as f:
            counts = load_is_counts(f)

    total = sum(counts.values())
    if total == 0:
        print("No IS data found.")
        return 1

    pcts = [float(p.strip()) for p in args.pcts.split(",") if p.strip()]
    pct_vals = percentiles(counts, pcts)

    arr, smoothed = smooth_counts(counts, args.smooth_window)
    peaks = detect_peaks(smoothed, args.peak_threshold)

    tail_count = sum(c for s, c in counts.items() if s >= args.tail)
    tail_pct = (tail_count / total) * 100.0

    mode_size = max(counts, key=lambda k: counts[k])
    mean = sum(size * count for size, count in counts.items()) / total

    print(f"Total pairs: {total}")
    print(f"Mode: {mode_size} bp")
    print(f"Mean: {mean:.2f} bp")
    for p in sorted(pcts):
        print(f"p{int(p)}: {pct_vals[p]} bp")
    print(f"Peaks (smoothed, threshold {args.peak_threshold:.2f}): {peaks}")
    print(f"Tail >= {args.tail} bp: {tail_count} ({tail_pct:.2f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
