#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  run_iceman_y_compare.sh [options]

Runs marker-state extraction + path ranking for multiple Iceman chrY callsets.
Intended to replace ad-hoc manual command sequences.

Options:
  --markers FILE        Marker VCF (default: auto-detect)
  --iceman-vcf FILE     DeepVariant VCF input (default: auto-detect)
  --iceman-gvcf FILE    DeepVariant gVCF input (default: auto-detect)
  --iceman-vcf-sample S Sample name for iceman.vcf (default: auto/first)
  --iceman-gvcf-sample S
                        Sample name for iceman.gvcf (default: auto/first)
  --solid-vcf FILE      Legacy SOLiD/freebayes VCF input
                        (default: auto-detect)
  --solid-sample S      Sample name for SOLiD VCF (default: Iceman_tst)
  --out-dir DIR         Output directory (default: /tmp/iceman_y_compare_YYYYmmdd_HHMMSS)
  --min-gq N            GQ threshold for derived.tsv rows (default: 20)
  --min-dp N            DP threshold for derived.tsv rows (default: 3)
  --skip-vcf            Skip iceman.vcf branch
  --skip-gvcf           Skip iceman.gvcf branch
  --skip-solid          Skip SOLiD branch
  -h, --help            Show this help

Outputs:
  <out-dir>/<dataset>.marker_status.tsv
  <out-dir>/<dataset>.derived.tsv
  <out-dir>/<dataset>.summary.txt
  <out-dir>/<dataset>.path_rank.tsv
  <out-dir>/key_markers.tsv
  <out-dir>/subtree_status.tsv
  <out-dir>/top_paths.tsv
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing required command: $1" >&2
    exit 1
  }
}

MARKERS="${YCOMPARE_MARKERS:-}"
ICEMAN_VCF="${YCOMPARE_ICEMAN_VCF:-}"
ICEMAN_GVCF="${YCOMPARE_ICEMAN_GVCF:-}"
ICEMAN_VCF_SAMPLE=""
ICEMAN_GVCF_SAMPLE=""
SOLID_VCF="${YCOMPARE_SOLID_VCF:-}"
SOLID_SAMPLE="Iceman_tst"
OUT_DIR="/tmp/iceman_y_compare_$(date +%Y%m%d_%H%M%S)"
MIN_GQ="20"
MIN_DP="3"
DO_VCF=1
DO_GVCF=1
DO_SOLID=1

while [ $# -gt 0 ]; do
  case "$1" in
    --markers) MARKERS="$2"; shift 2 ;;
    --iceman-vcf) ICEMAN_VCF="$2"; shift 2 ;;
    --iceman-gvcf) ICEMAN_GVCF="$2"; shift 2 ;;
    --iceman-vcf-sample) ICEMAN_VCF_SAMPLE="$2"; shift 2 ;;
    --iceman-gvcf-sample) ICEMAN_GVCF_SAMPLE="$2"; shift 2 ;;
    --solid-vcf) SOLID_VCF="$2"; shift 2 ;;
    --solid-sample) SOLID_SAMPLE="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --min-gq) MIN_GQ="$2"; shift 2 ;;
    --min-dp) MIN_DP="$2"; shift 2 ;;
    --skip-vcf) DO_VCF=0; shift ;;
    --skip-gvcf) DO_GVCF=0; shift ;;
    --skip-solid) DO_SOLID=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

need_cmd python3
need_cmd awk
need_cmd head

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

if [ -z "$MARKERS" ]; then
  for p in \
    "${PWD}/resources/snps_hg38.vcf.gz" \
    "${SCRIPT_DIR}/data/snps_hg38.vcf.gz" \
    "${PWD}/snps_hg38.vcf.gz"
  do
    if [ -f "$p" ]; then
      MARKERS="$p"
      break
    fi
  done
fi
[ -f "$MARKERS" ] || { echo "ERROR: markers not found. Set --markers FILE" >&2; exit 1; }

if [ -z "$ICEMAN_VCF" ]; then
  for p in "${PWD}/iceman.vcf" "${HOME}/iceman.vcf"; do
    if [ -f "$p" ]; then
      ICEMAN_VCF="$p"
      break
    fi
  done
fi
if [ -z "$ICEMAN_GVCF" ]; then
  for p in "${PWD}/iceman.gvcf" "${HOME}/iceman.gvcf"; do
    if [ -f "$p" ]; then
      ICEMAN_GVCF="$p"
      break
    fi
  done
fi
if [ -z "$SOLID_VCF" ]; then
  for p in \
    "${PWD}/chrY_called_Iceman_tst_hg38.vcf.gz"
  do
    if [ -f "$p" ]; then
      SOLID_VCF="$p"
      break
    fi
  done
fi

mkdir -p "$OUT_DIR"

run_one() {
  label="$1"
  input_vcf="$2"
  mode="$3"
  sample_name="$4"

  if [ ! -f "$input_vcf" ]; then
    echo "WARN: skipping $label, input missing: $input_vcf" >&2
    return 0
  fi

  out_prefix="${OUT_DIR}/${label}"
  echo "### [$label] input=$input_vcf mode=$mode"

  bgzip_flag=""
  if [ ! -f "${input_vcf}.tbi" ] && [ ! -f "${input_vcf}.csi" ]; then
    bgzip_flag="--bgzip-index-input"
  fi
  if [ -n "$sample_name" ]; then
    python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
      -i "$input_vcf" \
      --markers "$MARKERS" \
      -o "$out_prefix" \
      --sample "$sample_name" \
      --site-filter-mode "$mode" \
      --min-gq "$MIN_GQ" \
      --min-dp "$MIN_DP" \
      $bgzip_flag
  else
    python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
      -i "$input_vcf" \
      --markers "$MARKERS" \
      -o "$out_prefix" \
      --site-filter-mode "$mode" \
      --min-gq "$MIN_GQ" \
      --min-dp "$MIN_DP" \
      $bgzip_flag
  fi

  python3 "${SCRIPT_DIR}/y_path_rank.py" \
    --input "${out_prefix}.marker_status.tsv" \
    --out "${out_prefix}.path_rank.tsv" \
    --clade-prefix G \
    --normal-derived-weight 1 \
    --deam-derived-weight 0.5 \
    --ancestral-weight -1 \
    --nocall-weight 0
}

if [ "$DO_VCF" -eq 1 ]; then
  run_one "iceman_vcf" "$ICEMAN_VCF" "deepvariant" "$ICEMAN_VCF_SAMPLE"
fi
if [ "$DO_GVCF" -eq 1 ]; then
  run_one "iceman_gvcf" "$ICEMAN_GVCF" "deepvariant" "$ICEMAN_GVCF_SAMPLE"
fi
if [ "$DO_SOLID" -eq 1 ]; then
  run_one "iceman_solid" "$SOLID_VCF" "any" "$SOLID_SAMPLE"
fi

python3 - "$OUT_DIR" <<'PY'
import csv
import pathlib
import sys

out_dir = pathlib.Path(sys.argv[1])
datasets = []
for name in ("iceman_vcf", "iceman_gvcf", "iceman_solid"):
    p = out_dir / f"{name}.marker_status.tsv"
    if p.exists():
        datasets.append((name, p))

key_ids = {"M201", "L91", "L166", "Z6208", "Z31438", "Z6181", "Z6207", "Z6211"}
subtree_pat = "G2a2a1a2a1a1"

key_rows = []
subtree_rows = []
top_path_rows = []

for name, marker_file in datasets:
    status_counts = {}
    with marker_file.open("r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 13:
                continue
            marker_ids = cols[2].split(",")
            status = cols[12]
            isogg = cols[7]
            if any(m in key_ids for m in marker_ids):
                key_rows.append([
                    name, cols[2], cols[0], cols[1], status, cols[8], cols[9], cols[10], cols[7], cols[6]
                ])
            if subtree_pat in isogg:
                status_counts[status] = status_counts.get(status, 0) + 1

    for status, n in sorted(status_counts.items()):
        subtree_rows.append([name, subtree_pat, status, str(n)])

    path_file = out_dir / f"{name}.path_rank.tsv"
    if path_file.exists():
        with path_file.open("r", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            for i, row in enumerate(reader, start=1):
                if i > 10:
                    break
                top_path_rows.append([name] + row)

with (out_dir / "key_markers.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["dataset", "marker_id", "chrom", "pos", "status", "gt", "gq", "dp", "isogg", "hg"])
    w.writerows(sorted(key_rows, key=lambda r: (r[0], r[1], int(r[3]))))

with (out_dir / "subtree_status.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["dataset", "subtree", "status", "count"])
    w.writerows(subtree_rows)

with (out_dir / "top_paths.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([
        "dataset", "candidate", "score", "derived_total", "derived_transversion",
        "derived_deamination", "ancestral_total", "nocall_total", "ambiguous_total", "other_total"
    ])
    w.writerows(top_path_rows)
PY

echo "Done. Outputs in: $OUT_DIR"
