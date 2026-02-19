#!/usr/bin/env python3
"""
Marker-based Y haplogroup helper for both VCF and gVCF inputs.

This script streams marker rows (e.g. YBrowse/ISOGG) against sample calls and
classifies each marker as derived/ancestral/nocall/ambiguous.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


def is_gzip(path: Path) -> bool:
    with path.open("rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_text(path: Path):
    if is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def has_vcf_index(path: Path) -> bool:
    return path.with_suffix(path.suffix + ".tbi").exists() or path.with_suffix(path.suffix + ".csi").exists()


def line_stream_for_chrom(path: Path, chrom: str) -> Iterator[str]:
    """
    Stream VCF lines. Prefer indexed regional fetch via bcftools when available;
    fallback to full-file scan otherwise.
    """
    if shutil.which("bcftools") and has_vcf_index(path):
        proc = subprocess.Popen(
            ["bcftools", "view", "-r", chrom, "-Ov", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert proc.stdout is not None
        try:
            for ln in proc.stdout:
                yield ln
        finally:
            _out, err = proc.communicate()
            if proc.returncode not in (0, None):
                raise SystemExit(f"ERROR: bcftools regional fetch failed for {path}: {err.strip()}")
    else:
        with open_text(path) as fh:
            for ln in fh:
                yield ln


def parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not info or info == ".":
        return out
    for field in info.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = "1"
    return out


def site_filter_ok(filt: str, mode: str) -> bool:
    if mode == "any":
        return True
    if mode == "pass":
        return filt == "PASS"
    if mode == "pass-or-dot":
        return filt in {"PASS", "."}
    if mode == "deepvariant":
        return filt in {"PASS", "RefCall"}
    raise ValueError(f"unknown filter mode: {mode}")


@dataclass
class MarkerRow:
    chrom: str
    pos: int
    vid: str
    ref: str
    alt: str
    aa: str
    hg: str
    isogg: str


@dataclass
class SampleRow:
    chrom: str
    pos: int
    end: int
    filt: str
    ref: str
    alt: str
    gt: str
    gq: str
    dp: str
    is_block: bool


def parse_gt(gt: str) -> Tuple[str, Optional[int]]:
    if gt in {"", ".", "./.", ".|."}:
        return "nocall", None
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    alleles = [gt] if sep is None else gt.split(sep)
    if any(a == "." for a in alleles):
        return "nocall", None
    uniq = set(alleles)
    if len(uniq) != 1:
        return "ambiguous", None
    try:
        return "ok", int(next(iter(uniq)))
    except ValueError:
        return "ambiguous", None


def classify_marker(marker: MarkerRow, sample: Optional[SampleRow]) -> Tuple[str, str]:
    # returns (source, status)
    if sample is None:
        return "missing", "nocall"
    source = "gvcf_block" if sample.is_block else "variant"
    state, allele_idx = parse_gt(sample.gt)
    if state == "nocall":
        return source, "nocall"
    if state == "ambiguous":
        return source, "ambiguous"
    assert allele_idx is not None

    marker_ref = marker.ref
    marker_alt = marker.alt.split(",")[0]
    aa = marker.aa
    if not aa or aa == ".":
        return source, "unknown"

    # same logic as bcftools expression previously used:
    # derived if (GT=1 and AA=REF) or (GT=0 and AA=ALT)
    if allele_idx == 1 and aa == marker_ref:
        return source, "derived"
    if allele_idx == 0 and aa == marker_alt:
        return source, "derived"
    if allele_idx == 0 and aa == marker_ref:
        return source, "ancestral"
    if allele_idx == 1 and aa == marker_alt:
        return source, "ancestral"
    return source, "other"


def marker_rows(path: Path, chrom: str = "chrY") -> Iterator[MarkerRow]:
    for line in line_stream_for_chrom(path, chrom):
        if not line or line[0] == "#":
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8:
            continue
        if cols[0] != chrom:
            continue
        info = parse_info(cols[7])
        yield MarkerRow(
            chrom=cols[0],
            pos=int(cols[1]),
            vid=cols[2],
            ref=cols[3],
            alt=cols[4],
            aa=info.get("AA", "."),
            hg=info.get("HG", "."),
            isogg=info.get("ISOGG", "."),
        )


def sample_rows(path: Path, sample_name: Optional[str], chrom: str = "chrY") -> Iterator[SampleRow]:
    sample_idx: Optional[int] = None
    fmt_idx = 8
    for line in line_stream_for_chrom(path, chrom):
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t")
            if len(hdr) > 9:
                samples = hdr[9:]
                if sample_name:
                    if sample_name not in samples:
                        raise SystemExit(f"ERROR: sample '{sample_name}' not found in {path}")
                    sample_idx = 9 + samples.index(sample_name)
                else:
                    sample_idx = 9
            else:
                sample_idx = None
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8 or cols[0] != chrom:
            continue
        if sample_idx is None or len(cols) <= sample_idx:
            continue
        info = parse_info(cols[7])
        end = int(info.get("END", cols[1]))
        fmt_keys = cols[fmt_idx].split(":") if len(cols) > fmt_idx else []
        sample_vals = cols[sample_idx].split(":")
        sm = dict(zip(fmt_keys, sample_vals))
        alt = cols[4]
        is_block = "END" in info and ("<*>" in alt or alt == "<*>")
        yield SampleRow(
            chrom=cols[0],
            pos=int(cols[1]),
            end=end,
            filt=cols[6],
            ref=cols[3],
            alt=alt,
            gt=sm.get("GT", "."),
            gq=sm.get("GQ", "."),
            dp=sm.get("DP", sm.get("MIN_DP", ".")),
            is_block=is_block,
        )


def pick_variant(rows: List[SampleRow]) -> Optional[SampleRow]:
    if not rows:
        return None
    # Prefer non-block rows. For ties, prefer non-missing GT and higher GQ.
    non_block = [r for r in rows if not r.is_block]
    cands = non_block or rows

    def rank(row: SampleRow):
        gt_state, _ = parse_gt(row.gt)
        nonmiss = 1 if gt_state != "nocall" else 0
        try:
            gq = int(row.gq)
        except Exception:
            gq = -1
        return (nonmiss, gq)

    return sorted(cands, key=rank, reverse=True)[0]


def to_int_or_none(x: str) -> Optional[int]:
    try:
        return int(x)
    except Exception:
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Y marker status from VCF/gVCF.")
    ap.add_argument("-i", "--input", required=True, help="sample VCF/gVCF (plain or gzip)")
    ap.add_argument("--markers", required=True, help="marker VCF (YBrowse/ISOGG)")
    ap.add_argument("-o", "--out-prefix", required=True, help="output prefix")
    ap.add_argument("-s", "--sample", default="", help="sample name (default: first sample)")
    ap.add_argument("--chrom", default="chrY", help="target chromosome (default: chrY)")
    ap.add_argument("--site-filter-mode", default="any", choices=["any", "pass", "pass-or-dot", "deepvariant"])
    ap.add_argument("--min-gq", type=int, default=None)
    ap.add_argument("--min-dp", type=int, default=None)
    ap.add_argument(
        "--bgzip-index-input",
        action="store_true",
        help="When input has no index, create temporary bgzip+tabix copy for fast regional fetch.",
    )
    args = ap.parse_args()

    in_path = Path(args.input)
    marker_path = Path(args.markers)
    if not in_path.exists():
        raise SystemExit(f"ERROR: input not found: {in_path}")
    if not marker_path.exists():
        raise SystemExit(f"ERROR: markers not found: {marker_path}")

    temp_dir: Optional[tempfile.TemporaryDirectory[str]] = None
    sample_path_for_read = in_path
    if args.bgzip_index_input and not has_vcf_index(in_path):
        if not shutil.which("bgzip") or not shutil.which("tabix"):
            raise SystemExit("ERROR: --bgzip-index-input requires bgzip and tabix in PATH")
        temp_dir = tempfile.TemporaryDirectory(prefix="yhaplo_")
        sample_path_for_read = Path(temp_dir.name) / "sample.vcf.gz"
        # stream-convert (handles plain and gzip input)
        with sample_path_for_read.open("wb") as out_fh:
            p1 = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_fh)
            assert p1.stdin is not None
            src_fh = gzip.open(in_path, "rb") if is_gzip(in_path) else in_path.open("rb")
            try:
                with src_fh:
                    shutil.copyfileobj(src_fh, p1.stdin)
            finally:
                p1.stdin.close()
            rc = p1.wait()
        if rc != 0:
            raise SystemExit("ERROR: bgzip conversion failed")
        p2 = subprocess.Popen(
            ["tabix", "-f", "-p", "vcf", str(sample_path_for_read)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        _o2, e2 = p2.communicate()
        if p2.returncode != 0:
            raise SystemExit(f"ERROR: tabix indexing failed: {e2.strip()}")

    out_prefix = args.out_prefix
    marker_status_tsv = Path(f"{out_prefix}.marker_status.tsv")
    derived_tsv = Path(f"{out_prefix}.derived.tsv")
    summary_txt = Path(f"{out_prefix}.summary.txt")

    m_iter = marker_rows(marker_path, chrom=args.chrom)
    s_iter = sample_rows(sample_path_for_read, sample_name=(args.sample or None), chrom=args.chrom)

    cur_s: Optional[SampleRow] = next(s_iter, None)
    active_block: Optional[SampleRow] = None

    counts = {
        "total_markers": 0,
        "resolved_nonmiss": 0,
        "derived": 0,
        "ancestral": 0,
        "nocall": 0,
        "ambiguous": 0,
        "other": 0,
    }
    isogg_counts: Dict[str, int] = {}
    deepest: Dict[str, int] = {}
    max_depth = 0

    with marker_status_tsv.open("w", encoding="utf-8") as m_out, derived_tsv.open("w", encoding="utf-8") as d_out:
        for m in m_iter:
            counts["total_markers"] += 1
            pos = m.pos

            if active_block is not None and active_block.end < pos:
                active_block = None

            while cur_s is not None and cur_s.pos < pos:
                if site_filter_ok(cur_s.filt, args.site_filter_mode):
                    if cur_s.is_block and cur_s.end >= pos:
                        active_block = cur_s
                cur_s = next(s_iter, None)

            same_pos: List[SampleRow] = []
            while cur_s is not None and cur_s.pos == pos:
                if site_filter_ok(cur_s.filt, args.site_filter_mode):
                    same_pos.append(cur_s)
                    if cur_s.is_block and cur_s.end >= pos:
                        active_block = cur_s
                cur_s = next(s_iter, None)

            picked = pick_variant(same_pos)
            if picked is None and active_block is not None and active_block.end >= pos and site_filter_ok(
                active_block.filt, args.site_filter_mode
            ):
                picked = active_block

            source, status = classify_marker(m, picked)

            gt = picked.gt if picked else "."
            gq = picked.gq if picked else "."
            dp = picked.dp if picked else "."

            m_out.write(
                "\t".join(
                    [
                        m.chrom,
                        str(m.pos),
                        m.vid,
                        m.ref,
                        m.alt,
                        m.aa,
                        m.hg,
                        m.isogg,
                        gt,
                        gq,
                        dp,
                        source,
                        status,
                    ]
                )
                + "\n"
            )

            if status != "nocall":
                counts["resolved_nonmiss"] += 1
            counts[status] = counts.get(status, 0) + 1

            if status == "derived":
                gq_i = to_int_or_none(gq)
                dp_i = to_int_or_none(dp)
                if args.min_gq is not None and (gq_i is None or gq_i < args.min_gq):
                    continue
                if args.min_dp is not None and (dp_i is None or dp_i < args.min_dp):
                    continue
                d_out.write(
                    "\t".join(
                        [
                            m.chrom,
                            str(m.pos),
                            m.vid,
                            m.ref,
                            m.alt,
                            m.aa,
                            m.hg,
                            m.isogg,
                            gt,
                            gq,
                            dp,
                            source,
                            status,
                        ]
                    )
                    + "\n"
                )
                counts["derived"] += 0  # already counted above; keep no-op for clarity

                isogg = m.isogg.strip('"') if m.isogg and m.isogg != "." else "(missing)"
                isogg_counts[isogg] = isogg_counts.get(isogg, 0) + 1
                hg = m.hg.strip('"')
                if hg and hg != ".":
                    d = hg.count(">") + 1
                    if d > max_depth:
                        max_depth = d
                        deepest = {hg: 1}
                    elif d == max_depth:
                        deepest[hg] = 1

    with summary_txt.open("w", encoding="utf-8") as out:
        out.write(f"Input: {in_path}\n")
        out.write(f"Markers: {marker_path}\n")
        out.write(f"Chrom: {args.chrom}\n")
        out.write(f"FilterMode: {args.site_filter_mode}\n")
        out.write(f"MinGQ: {args.min_gq if args.min_gq is not None else 'none'}\n")
        out.write(f"MinDP: {args.min_dp if args.min_dp is not None else 'none'}\n\n")
        out.write("Counts:\n")
        for k in ["total_markers", "resolved_nonmiss", "derived", "ancestral", "nocall", "ambiguous", "other"]:
            out.write(f"  {k}: {counts.get(k, 0)}\n")
        out.write("\nDeepest HG candidates (derived rows only):\n")
        if max_depth == 0:
            out.write("  (none)\n")
        else:
            out.write(f"  depth={max_depth}\n")
            for h in sorted(deepest):
                out.write(f"  {h}\n")
        out.write("\nTop ISOGG labels among derived rows:\n")
        for label, n in sorted(isogg_counts.items(), key=lambda kv: (-kv[1], kv[0]))[:30]:
            out.write(f"  {n}\t{label}\n")

    if temp_dir is not None:
        temp_dir.cleanup()
    print(f"Done. marker_status={marker_status_tsv} derived={derived_tsv} summary={summary_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
