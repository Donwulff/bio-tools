#!/bin/sh
set -eu

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
    [--java-opts "-Xms1g -Xmx8g"] \
    [--out-dir DIR]

Purpose:
  Runs both Y-haplogroup paths for cross-check:
  1) hg38 markers -> hs1 liftover, then call against native hs1 sample
  2) hs1 sample VCF -> hg38 liftover, then call against hg38 markers
  3) (if available) rank hg38-lifted sample against hg38 marker GFF3 for YFull-style terminal labels
  4) rank candidate paths (auto top-level clade by default)

Defaults:
  out-dir defaults to ./experiments/y_dual_liftover_<timestamp>
  markers default is fetched to ./resources/snps_hg38.vcf.gz
  refs auto-detect from --ref-root, BIO_TOOLS_INDEX_DIR/INDEX_DIR, ./resources/ref, mapping/index

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
MARKERS_HG38=""
MARKERS_GFF3_HG38=""
PICARD_JAR=""
OUT_DIR=""
JAVA_OPTS=""
ALLOW_EXTENDED_REF=0
CLADE_PREFIX="${YPATH_CLADE_PREFIX:-R}"
AUTO_CLADE=1

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
    --java-opts) JAVA_OPTS="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
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

need_cmd java
need_cmd bcftools
need_cmd tabix
need_cmd python3
need_cmd samtools

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "${SCRIPT_DIR}/lib_liftover.sh"

# Optional local config (same convention as mapping scripts).
if [ -f "${BIO_TOOLS_CFG:-./bio-tools.cfg}" ]; then
  # shellcheck disable=SC1090
  . "${BIO_TOOLS_CFG:-./bio-tools.cfg}"
fi

if [ -z "$REF_ROOT" ]; then
  REF_ROOT="${BIO_TOOLS_INDEX_DIR:-${INDEX_DIR:-}}"
fi

pick_ref() {
  roots="$1"
  shift
  for r in $roots; do
    [ -n "$r" ] || continue
    for n in "$@"; do
      p="${r%/}/$n"
      if [ -f "$p" ]; then
        echo "$p"
        return 0
      fi
    done
  done
  return 1
}

if [ -z "$REF_HS1" ] || [ -z "$REF_HG38" ]; then
  roots=""
  [ -n "$REF_ROOT" ] && roots="$roots $REF_ROOT"
  roots="$roots ${PWD}/resources/ref ${SCRIPT_DIR}/../mapping/index"

  if [ -z "$REF_HS1" ]; then
    REF_HS1="$(pick_ref "$roots" \
      "chm13v2.0_maskedY_rCRS.fa.gz" \
      "chm13v2.0_maskedY_rCRS.H3630.fa.gz" \
      "chm13v2.0_maskedY_rCRSDH3630O1102.fa" \
      "chm13v2.0_maskedY_rCRS.fa" \
      "chm13v2.0_maskedY_rCRS.H3630.fa" || true)"
  fi
  if [ -z "$REF_HG38" ]; then
    REF_HG38="$(pick_ref "$roots" \
      "GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna" \
      "hg38p14.H3630.fa.gz" \
      "hg38p14DH3630O1102.fa" \
      "hg38p14.H3630.fa" || true)"
  fi
fi

[ -n "$REF_HS1" ] || { echo "ERROR: --ref-hs1 not set and auto-detect failed." >&2; exit 1; }
[ -n "$REF_HG38" ] || { echo "ERROR: --ref-hg38 not set and auto-detect failed." >&2; exit 1; }
[ -f "$REF_HS1" ] || { echo "ERROR: not found: $REF_HS1" >&2; exit 1; }
[ -f "$REF_HG38" ] || { echo "ERROR: not found: $REF_HG38" >&2; exit 1; }

if [ "$ALLOW_EXTENDED_REF" -eq 0 ]; then
  case "$(basename "$REF_HS1") $(basename "$REF_HG38")" in
    *H3630*|*DH3630*|*O1102*|*hs38d1*|*oral_microbiome*)
      echo "ERROR: selected reference appears to include extended/decoy content (H3630/hs38d1/O1102)." >&2
      echo "This workflow defaults to non-extended references for Y analysis." >&2
      echo "Pick non-extended refs explicitly, or pass --allow-extended-ref to override." >&2
      exit 2
      ;;
  esac
fi

echo "Using refs:"
echo "  hs1  $REF_HS1"
echo "  hg38 $REF_HG38"

if [ -z "$PICARD_JAR" ]; then
  PICARD_JAR="$(detect_picard_jar "$SCRIPT_DIR" || true)"
fi
[ -n "$PICARD_JAR" ] || { echo "ERROR: set --picard-jar (picard.jar not auto-detected)" >&2; exit 1; }
[ -f "$PICARD_JAR" ] || { echo "ERROR: not found: $PICARD_JAR" >&2; exit 1; }
if [ -z "$JAVA_OPTS" ]; then
  JAVA_OPTS="${PICARD_JAVA_OPTS:-"-Xms1g -Xmx8g"}"
fi

ensure_ref_sidecars "$REF_HS1" "$PICARD_JAR" "$JAVA_OPTS"
ensure_ref_sidecars "$REF_HG38" "$PICARD_JAR" "$JAVA_OPTS"

if [ -z "$MARKERS_HG38" ]; then
  if [ -n "$MARKERS_GFF3_HG38" ]; then
    sib="${MARKERS_GFF3_HG38%.*}.vcf.gz"
    if [ -f "$sib" ]; then
      MARKERS_HG38="$sib"
    fi
  fi
fi
if [ -z "$MARKERS_HG38" ]; then
  if [ -n "${Y_RESOURCES_DIR:-}" ] && [ -f "${Y_RESOURCES_DIR%/}/snps_hg38.vcf.gz" ]; then
    MARKERS_HG38="${Y_RESOURCES_DIR%/}/snps_hg38.vcf.gz"
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
java $JAVA_OPTS -jar "$PICARD_JAR" LiftoverVcf \
  I="$MARKERS_HG38_CLEAN" \
  O="$MARKERS_HS1" \
  CHAIN="$CHAIN_HG38_TO_HS1" \
  R="$REF_HS1" \
  REJECT="$MARKERS_HS1_REJECT" \
  RECOVER_SWAPPED_REF_ALT=true \
  WARN_ON_MISSING_CONTIG=true
tabix -f "$MARKERS_HS1"

echo "Liftover sample VCF: hs1 -> hg38"
java $JAVA_OPTS -jar "$PICARD_JAR" LiftoverVcf \
  I="$SAMPLE_VCF_HS1" \
  O="$SAMPLE_HG38" \
  CHAIN="$CHAIN_HS1_TO_HG38" \
  R="$REF_HG38" \
  REJECT="$SAMPLE_HG38_REJECT" \
  RECOVER_SWAPPED_REF_ALT=true \
  WARN_ON_MISSING_CONTIG=true
tabix -f "$SAMPLE_HG38"

echo "Run marker-state extraction on hs1 branch"
python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
  -i "$SAMPLE_VCF_HS1" \
  --markers "$MARKERS_HS1" \
  -o "${OUT_DIR}/hs1_from_markers_lifted" \
  --site-filter-mode deepvariant \
  --bgzip-index-input

if [ -n "$SAMPLE_GVCF_HS1" ]; then
  python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
    -i "$SAMPLE_GVCF_HS1" \
    --markers "$MARKERS_HS1" \
    -o "${OUT_DIR}/hs1_from_markers_lifted_gvcf" \
    --site-filter-mode deepvariant \
    --bgzip-index-input
fi

echo "Run marker-state extraction on hg38-lifted sample branch"
python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
  -i "$SAMPLE_HG38" \
  --markers "$MARKERS_HG38_CLEAN" \
  -o "${OUT_DIR}/hg38_from_sample_lifted" \
  --site-filter-mode deepvariant \
  --bgzip-index-input

if [ -n "$MARKERS_GFF3_HG38" ] && [ -f "$MARKERS_GFF3_HG38" ]; then
  echo "Run marker-state extraction on hg38-lifted sample branch (GFF3 markers)"
  python3 "${SCRIPT_DIR}/y_haplo_from_markers.py" \
    -i "$SAMPLE_HG38" \
    --markers "$MARKERS_GFF3_HG38" \
    -o "${OUT_DIR}/hg38_from_sample_lifted_gff3" \
    --site-filter-mode deepvariant \
    --bgzip-index-input
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
