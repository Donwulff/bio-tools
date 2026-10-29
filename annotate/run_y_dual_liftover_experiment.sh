#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
. "${SCRIPT_DIR}/lib_liftover.sh"
load_liftover_config

usage() {
  cat <<'EOF'
Usage:
  run_y_dual_liftover_experiment.sh \
    --sample-vcf-hs1 SAMPLE_HS1.vcf.gz \
    [--sample-gvcf-hs1 SAMPLE_HS1.gvcf.gz] \
    [--ref-hs1 REF_HS1.fa] \
    [--ref-hg38 REF_HG38.fa] \
    [--ref-root DIR] \
    --chain-hg38-to-hs1 hg38ToHs1.chain.gz \
    --chain-hs1-to-hg38 hs1ToHg38.chain.gz \
    [--auto-clade] [--no-auto-clade] [--clade-prefix PREFIX] \
    [--markers-hg38 snps_hg38.vcf.gz] \
    [--markers-gff3-hg38 snps_hg38.gff3] \
    [--allow-extended-ref] \
    [--picard-jar picard.jar] \
    [--backend auto|gatk|picard-jar|picard-cli] \
    [--java-opts "-Xms1g -Xmx8g"] \
    [--min-gq N] [--min-dp N] \
    [--out-dir DIR]

Purpose:
  Runs both Y-haplogroup paths for cross-check:
  1) hg38 markers -> hs1 liftover, then call against native hs1 sample
  2) hs1 sample VCF -> hg38 liftover, then call against hg38 markers
  3) (if available) rank hg38-lifted sample against hg38 marker GFF3 for YFull-style terminal labels
  4) rank candidate paths (auto top-level clade by default)

Defaults:
  out-dir defaults to ./experiments/y_dual_liftover_<timestamp>
  markers default to Y_RESOURCES_DIR or ./resources (explicit paths win)
  refs share discovery with the batch driver; non-extended references are preferred
  trusted BIO_TOOLS_CFG or ./bio-tools.cfg loads before defaults and CLI options
  quality thresholds default to unset; when supplied they apply to all marker states
  reference-based marker resolution requires pysam

Output:
  <out-dir>/summary.tsv
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

SAMPLE_VCF_HS1=""
SAMPLE_GVCF_HS1=""
REF_HS1="${REF_HS1:-}"
REF_HG38="${REF_HG38:-}"
REF_ROOT="${REF_ROOT:-}"
CHAIN_HG38_TO_HS1=""
CHAIN_HS1_TO_HG38=""
MARKERS_HG38="${YMARKERS_OUT:-}"
MARKERS_GFF3_HG38="${YMARKERS_GFF3_OUT:-}"
PICARD_JAR="${PICARD_JAR:-}"
LIFTOVER_BACKEND="${LIFTOVER_BACKEND:-auto}"
OUT_DIR=""
LIFTOVER_JAVA_OPTS="${PICARD_JAVA_OPTS-"-Xms1g -Xmx8g"}"
ALLOW_EXTENDED_REF=0
CLADE_PREFIX="${YPATH_CLADE_PREFIX:-R}"
AUTO_CLADE=1
MIN_GQ=""
MIN_DP=""

while [ $# -gt 0 ]; do
  case "$1" in
    --sample-vcf-hs1) SAMPLE_VCF_HS1="$2"; shift 2 ;;
    --sample-gvcf-hs1) SAMPLE_GVCF_HS1="$2"; shift 2 ;;
    --ref-hs1) REF_HS1="$2"; shift 2 ;;
    --ref-hg38) REF_HG38="$2"; shift 2 ;;
    --ref-root) REF_ROOT="$2"; shift 2 ;;
    --chain-hg38-to-hs1) CHAIN_HG38_TO_HS1="$2"; shift 2 ;;
    --chain-hs1-to-hg38) CHAIN_HS1_TO_HG38="$2"; shift 2 ;;
    --auto-clade) AUTO_CLADE=1; shift ;;
    --no-auto-clade) AUTO_CLADE=0; shift ;;
    --clade-prefix) CLADE_PREFIX="$2"; AUTO_CLADE=0; shift 2 ;;
    --markers-hg38) MARKERS_HG38="$2"; shift 2 ;;
    --markers-gff3-hg38) MARKERS_GFF3_HG38="$2"; shift 2 ;;
    --allow-extended-ref) ALLOW_EXTENDED_REF=1; shift ;;
    --picard-jar) PICARD_JAR="$2"; shift 2 ;;
    --backend) LIFTOVER_BACKEND="$2"; shift 2 ;;
    --java-opts) LIFTOVER_JAVA_OPTS="$2"; shift 2 ;;
    --min-gq) MIN_GQ="$2"; shift 2 ;;
    --min-dp) MIN_DP="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

for threshold in "$MIN_GQ" "$MIN_DP"; do
  case "$threshold" in
    *[!0-9]*) echo "ERROR: --min-gq and --min-dp must be nonnegative integers" >&2; exit 1 ;;
  esac
done

[ -n "$SAMPLE_VCF_HS1" ] || { echo "ERROR: --sample-vcf-hs1 is required" >&2; exit 1; }
[ -n "$CHAIN_HG38_TO_HS1" ] || { echo "ERROR: --chain-hg38-to-hs1 is required" >&2; exit 1; }
[ -n "$CHAIN_HS1_TO_HG38" ] || { echo "ERROR: --chain-hs1-to-hg38 is required" >&2; exit 1; }

[ -f "$SAMPLE_VCF_HS1" ] || { echo "ERROR: not found: $SAMPLE_VCF_HS1" >&2; exit 1; }
[ -f "$CHAIN_HG38_TO_HS1" ] || { echo "ERROR: not found: $CHAIN_HG38_TO_HS1" >&2; exit 1; }
[ -f "$CHAIN_HS1_TO_HG38" ] || { echo "ERROR: not found: $CHAIN_HS1_TO_HG38" >&2; exit 1; }

if [ -n "$SAMPLE_GVCF_HS1" ]; then
  [ -f "$SAMPLE_GVCF_HS1" ] || { echo "ERROR: not found: $SAMPLE_GVCF_HS1" >&2; exit 1; }
fi

need_cmd bcftools
need_cmd tabix
need_cmd python3
need_cmd samtools
python3 -c 'import pysam' >/dev/null 2>&1 || {
  echo "ERROR: reference-based marker resolution requires pysam (python3 -m pip install pysam)" >&2
  exit 1
}

REF_HS1="$(find_liftover_ref hs1 "$REF_HS1" "$REF_ROOT" "$ALLOW_EXTENDED_REF")"
REF_HG38="$(find_liftover_ref hg38 "$REF_HG38" "$REF_ROOT" "$ALLOW_EXTENDED_REF")"

echo "Using refs:"
echo "  hs1  $REF_HS1"
echo "  hg38 $REF_HG38"

BACKEND_SPEC="$(select_liftover_backend "$LIFTOVER_BACKEND" "$PICARD_JAR" "$LIFTOVER_JAVA_OPTS")"

ensure_ref_sidecars "$REF_HS1" "$PICARD_JAR" "$LIFTOVER_JAVA_OPTS" "$BACKEND_SPEC"
ensure_ref_sidecars "$REF_HG38" "$PICARD_JAR" "$LIFTOVER_JAVA_OPTS" "$BACKEND_SPEC"

if [ -z "$MARKERS_HG38" ]; then
  if [ -n "$MARKERS_GFF3_HG38" ]; then
    sib="${MARKERS_GFF3_HG38%.*}.vcf.gz"
    if [ -f "$sib" ]; then
      MARKERS_HG38="$sib"
    fi
  fi
fi
if [ -z "$MARKERS_HG38" ]; then
  if [ -n "${Y_RESOURCES_DIR:-}" ]; then
    MARKERS_HG38="${Y_RESOURCES_DIR%/}/snps_hg38.vcf.gz"
    if [ ! -f "$MARKERS_HG38" ]; then
      if [ -n "$MARKERS_GFF3_HG38" ]; then
        "${SCRIPT_DIR}/fetch_ybrowse_markers.sh" --out "$MARKERS_HG38" --no-gff3
      else
        "${SCRIPT_DIR}/fetch_ybrowse_markers.sh" --out "$MARKERS_HG38"
      fi
    fi
  fi
fi
if [ -z "$MARKERS_HG38" ]; then
  for p in \
    "${PWD}/resources/snps_hg38.vcf.gz" \
    "${PWD}/snps_hg38.vcf.gz" \
    "${SCRIPT_DIR}/snps_hg38.vcf.gz"
  do
    if [ -f "$p" ]; then
      MARKERS_HG38="$p"
      break
    fi
  done
fi
if [ -z "$MARKERS_HG38" ]; then
  MARKERS_HG38="$(pwd)/resources/snps_hg38.vcf.gz"
  "${SCRIPT_DIR}/fetch_ybrowse_markers.sh" --out "$MARKERS_HG38"
fi
[ -f "$MARKERS_HG38" ] || { echo "ERROR: not found: $MARKERS_HG38" >&2; exit 1; }
[ -f "${MARKERS_HG38}.tbi" ] || tabix -f "$MARKERS_HG38"

if [ -z "$MARKERS_GFF3_HG38" ]; then
  if [ -n "${Y_RESOURCES_DIR:-}" ] && [ -f "${Y_RESOURCES_DIR%/}/snps_hg38.gff3" ]; then
    MARKERS_GFF3_HG38="${Y_RESOURCES_DIR%/}/snps_hg38.gff3"
  fi
fi
if [ -z "$MARKERS_GFF3_HG38" ]; then
  for p in \
    "${PWD}/resources/snps_hg38.gff3" \
    "${PWD}/snps_hg38.gff3" \
    "${SCRIPT_DIR}/snps_hg38.gff3"
  do
    if [ -f "$p" ]; then
      MARKERS_GFF3_HG38="$p"
      break
    fi
  done
fi

if [ -z "$OUT_DIR" ]; then
  OUT_DIR="$(pwd)/experiments/y_dual_liftover_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUT_DIR"

MARKERS_HS1="${OUT_DIR}/markers_hs1.vcf.gz"
MARKERS_HS1_REJECT="${OUT_DIR}/markers_hs1.reject.vcf.gz"
SAMPLE_HG38="${OUT_DIR}/sample_hg38_lifted.vcf.gz"
SAMPLE_HG38_REJECT="${OUT_DIR}/sample_hg38_lifted.reject.vcf.gz"
MARKERS_HG38_CLEAN="${OUT_DIR}/markers_hg38.clean.vcf.gz"

"${SCRIPT_DIR}/sanitize_marker_vcf.sh" --input "$MARKERS_HG38" --output "$MARKERS_HG38_CLEAN"

echo "Liftover markers: hg38 -> hs1"
run_liftover_tool "$BACKEND_SPEC" LiftoverVcf \
  --INPUT "$MARKERS_HG38_CLEAN" \
  --OUTPUT "$MARKERS_HS1" \
  --CHAIN "$CHAIN_HG38_TO_HS1" \
  --REFERENCE_SEQUENCE "$REF_HS1" \
  --REJECT "$MARKERS_HS1_REJECT" \
  --RECOVER_SWAPPED_REF_ALT true \
  --WARN_ON_MISSING_CONTIG true
tabix -f "$MARKERS_HS1"

echo "Liftover sample VCF: hs1 -> hg38"
run_liftover_tool "$BACKEND_SPEC" LiftoverVcf \
  --INPUT "$SAMPLE_VCF_HS1" \
  --OUTPUT "$SAMPLE_HG38" \
  --CHAIN "$CHAIN_HS1_TO_HG38" \
  --REFERENCE_SEQUENCE "$REF_HG38" \
  --REJECT "$SAMPLE_HG38_REJECT" \
  --RECOVER_SWAPPED_REF_ALT true \
  --WARN_ON_MISSING_CONTIG true
tabix -f "$SAMPLE_HG38"

echo "Run marker-state extraction on hs1 branch"
extract_markers() {
  set -- -i "$1" --markers "$2" -o "$3" --reference "$4" \
    --site-filter-mode deepvariant --bgzip-index-input
  [ -z "$MIN_GQ" ] || set -- "$@" --min-gq "$MIN_GQ"
  [ -z "$MIN_DP" ] || set -- "$@" --min-dp "$MIN_DP"
  python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" "$@"
}

extract_markers "$SAMPLE_VCF_HS1" "$MARKERS_HS1" \
  "${OUT_DIR}/hs1_from_markers_lifted" "$REF_HS1"

if [ -n "$SAMPLE_GVCF_HS1" ]; then
  extract_markers "$SAMPLE_GVCF_HS1" "$MARKERS_HS1" \
    "${OUT_DIR}/hs1_from_markers_lifted_gvcf" "$REF_HS1"
fi

echo "Run marker-state extraction on hg38-lifted sample branch"
extract_markers "$SAMPLE_HG38" "$MARKERS_HG38_CLEAN" \
  "${OUT_DIR}/hg38_from_sample_lifted" "$REF_HG38"

if [ -n "$MARKERS_GFF3_HG38" ] && [ -f "$MARKERS_GFF3_HG38" ]; then
  echo "Run marker-state extraction on hg38-lifted sample branch (GFF3 markers)"
  extract_markers "$SAMPLE_HG38" "$MARKERS_GFF3_HG38" \
    "${OUT_DIR}/hg38_from_sample_lifted_gff3" "$REF_HG38"
fi

if [ "$AUTO_CLADE" -eq 1 ]; then
  python3 "${SCRIPT_DIR}/y_path_rank.py" \
    --input "${OUT_DIR}/hs1_from_markers_lifted.marker_status.tsv" \
    --out "${OUT_DIR}/hs1.paths.tsv" \
    --auto-clade
  python3 "${SCRIPT_DIR}/y_path_rank.py" \
    --input "${OUT_DIR}/hg38_from_sample_lifted.marker_status.tsv" \
    --out "${OUT_DIR}/hg38.paths.tsv" \
    --auto-clade
  if [ -f "${OUT_DIR}/hg38_from_sample_lifted_gff3.marker_status.tsv" ]; then
    python3 "${SCRIPT_DIR}/y_path_rank.py" \
      --input "${OUT_DIR}/hg38_from_sample_lifted_gff3.marker_status.tsv" \
      --out "${OUT_DIR}/hg38_gff3.paths.tsv" \
      --auto-clade
  fi
else
  python3 "${SCRIPT_DIR}/y_path_rank.py" \
    --input "${OUT_DIR}/hs1_from_markers_lifted.marker_status.tsv" \
    --out "${OUT_DIR}/hs1.paths.tsv" \
    --clade-prefix "$CLADE_PREFIX"
  python3 "${SCRIPT_DIR}/y_path_rank.py" \
    --input "${OUT_DIR}/hg38_from_sample_lifted.marker_status.tsv" \
    --out "${OUT_DIR}/hg38.paths.tsv" \
    --clade-prefix "$CLADE_PREFIX"
  if [ -f "${OUT_DIR}/hg38_from_sample_lifted_gff3.marker_status.tsv" ]; then
    python3 "${SCRIPT_DIR}/y_path_rank.py" \
      --input "${OUT_DIR}/hg38_from_sample_lifted_gff3.marker_status.tsv" \
      --out "${OUT_DIR}/hg38_gff3.paths.tsv" \
      --clade-prefix "$CLADE_PREFIX"
  fi
fi

{
  printf '%s\n' "dataset	candidate	score	derived_total	derived_transversion	derived_deamination	ancestral_total	nocall_total	ambiguous_total	other_total	net_support	support_ratio	derived_dp_sum	derived_dp_mean	derived_dp_max	derived_gq_mean	label_style"
  for p in hs1.paths.tsv hg38.paths.tsv hg38_gff3.paths.tsv; do
    f="${OUT_DIR}/${p}"
    [ -f "$f" ] || continue
    ds="${p%.paths.tsv}"
    awk -v d="$ds" 'BEGIN{FS=OFS="\t"} NR==1{next} ($2+0)>0 {print d,$0}' "$f" \
      | sed -n '1,30p'
  done
} > "${OUT_DIR}/top_paths.tsv"

{
  printf 'branch\tmetric\tvalue\n'
  printf 'global\tauto_clade\t%s\n' "${AUTO_CLADE}"
  printf 'global\tclade_prefix\t%s\n' "${CLADE_PREFIX}"
  awk '
    BEGIN{b="hs1"} /^  derived:/{print b"\tderived\t"$2}
    /^  ancestral:/{print b"\tancestral\t"$2}
    /^  nocall:/{print b"\tnocall\t"$2}
  ' "${OUT_DIR}/hs1_from_markers_lifted.summary.txt"
  awk '
    BEGIN{b="hg38_lifted_sample"} /^  derived:/{print b"\tderived\t"$2}
    /^  ancestral:/{print b"\tancestral\t"$2}
    /^  nocall:/{print b"\tnocall\t"$2}
  ' "${OUT_DIR}/hg38_from_sample_lifted.summary.txt"
  if [ -f "${OUT_DIR}/hg38_from_sample_lifted_gff3.summary.txt" ]; then
    awk '
      BEGIN{b="hg38_lifted_sample_gff3"} /^  derived:/{print b"\tderived\t"$2}
      /^  ancestral:/{print b"\tancestral\t"$2}
      /^  nocall:/{print b"\tnocall\t"$2}
    ' "${OUT_DIR}/hg38_from_sample_lifted_gff3.summary.txt"
  fi
} > "${OUT_DIR}/summary.tsv"

echo "Done. Output: $OUT_DIR"
echo "Key file: ${OUT_DIR}/summary.tsv"
echo "Top ranked candidates: ${OUT_DIR}/top_paths.tsv"
