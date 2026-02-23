#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  run_modern_y_experiment.sh --sample-vcf sample.vcf.gz [--sample-gvcf sample.gvcf.gz]
                             [--markers markers.vcf.gz] [--out-dir DIR]
                             [--allow-assembly-mismatch]

Purpose:
  Minimal driver for modern-sample chrY haplogroup experiments.
  - Detects chrY assembly mismatch between sample and marker VCF
  - Runs marker-state extraction when compatible

Notes:
  - Default marker path prefers /mnt/GenomicData/Iceman_haplo/snps_hg38.vcf.gz,
    then annotate/data/snps_hg38.vcf.gz.
  - For CHM13 sample calls and hg38 markers, liftover markers first or use
    a native GRCh38 callset.
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

SAMPLE_VCF=""
SAMPLE_GVCF=""
MARKERS=""
OUT_DIR=""
ALLOW_MISMATCH=0

while [ $# -gt 0 ]; do
  case "$1" in
    --sample-vcf) SAMPLE_VCF="$2"; shift 2 ;;
    --sample-gvcf) SAMPLE_GVCF="$2"; shift 2 ;;
    --markers) MARKERS="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --allow-assembly-mismatch) ALLOW_MISMATCH=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$SAMPLE_VCF" ] || { echo "ERROR: --sample-vcf is required" >&2; usage; exit 1; }
[ -f "$SAMPLE_VCF" ] || { echo "ERROR: sample VCF not found: $SAMPLE_VCF" >&2; exit 1; }

need_cmd bcftools
need_cmd python3

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)

if [ -z "$MARKERS" ]; then
  if [ -f "/mnt/GenomicData/Iceman_haplo/snps_hg38.vcf.gz" ]; then
    MARKERS="/mnt/GenomicData/Iceman_haplo/snps_hg38.vcf.gz"
  else
    MARKERS="${SCRIPT_DIR}/data/snps_hg38.vcf.gz"
  fi
fi
[ -f "$MARKERS" ] || {
  echo "ERROR: marker VCF not found: $MARKERS" >&2
  echo "Provide --markers PATH." >&2
  exit 1
}

if [ -z "$OUT_DIR" ]; then
  base=$(basename "$SAMPLE_VCF")
  base=${base%.vcf.gz}
  base=${base%.vcf}
  OUT_DIR="${REPO_ROOT}/mapping/tst/modern_y_${base}"
fi
mkdir -p "$OUT_DIR"

get_chr_len() {
  vcf="$1"
  chrom="$2"
  bcftools view -h "$vcf" 2>/dev/null \
    | awk -v c="$chrom" '
      $0 ~ /^##contig=<ID=/ {
        id=""; len="";
        if (match($0, /ID=[^,>]+/)) id=substr($0, RSTART+3, RLENGTH-3);
        if (match($0, /length=[0-9]+/)) len=substr($0, RSTART+7, RLENGTH-7);
        if (id == c) { print len; exit }
      }
    '
}

sample_y_len=$(get_chr_len "$SAMPLE_VCF" "chrY")
marker_y_len=$(get_chr_len "$MARKERS" "chrY")

if [ -z "$sample_y_len" ] || [ -z "$marker_y_len" ]; then
  echo "ERROR: could not determine chrY contig length from headers." >&2
  echo "sample chrY length: ${sample_y_len:-missing}" >&2
  echo "marker chrY length: ${marker_y_len:-missing}" >&2
  exit 1
fi

{
  echo "sample_vcf=$SAMPLE_VCF"
  echo "sample_gvcf=${SAMPLE_GVCF:-none}"
  echo "markers=$MARKERS"
  echo "sample_chrY_length=$sample_y_len"
  echo "marker_chrY_length=$marker_y_len"
} > "${OUT_DIR}/run.info"

if [ "$sample_y_len" != "$marker_y_len" ] && [ "$ALLOW_MISMATCH" -eq 0 ]; then
  cat > "${OUT_DIR}/STOP_ASSEMBLY_MISMATCH.txt" <<EOF
Assembly mismatch detected (chrY lengths differ):
  sample chrY length: ${sample_y_len}
  marker chrY length: ${marker_y_len}

This run was stopped intentionally to avoid invalid marker assignment.

Next step:
  - liftover markers to sample assembly, OR
  - run haplogroup assignment on a matching-assembly callset.
EOF
  echo "Stopped due to assembly mismatch. See: ${OUT_DIR}/STOP_ASSEMBLY_MISMATCH.txt"
  exit 2
fi

if [ -n "$SAMPLE_GVCF" ]; then
  [ -f "$SAMPLE_GVCF" ] || { echo "ERROR: sample gVCF not found: $SAMPLE_GVCF" >&2; exit 1; }
  python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
    -i "$SAMPLE_GVCF" \
    --markers "$MARKERS" \
    -o "${OUT_DIR}/from_gvcf" \
    --site-filter-mode deepvariant \
    --bgzip-index-input
fi

python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
  -i "$SAMPLE_VCF" \
  --markers "$MARKERS" \
  -o "${OUT_DIR}/from_vcf" \
  --site-filter-mode deepvariant \
  --bgzip-index-input

echo "Done. Outputs under: $OUT_DIR"
