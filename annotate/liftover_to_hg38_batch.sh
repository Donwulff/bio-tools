#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  liftover_to_hg38_batch.sh \
    --chain-hs1-to-hg38 hs1ToHg38.chain.gz \
    [--ref-hg38 GRCh38.fa] \
    [--ref-root DIR] \
    [--picard-jar picard.jar] \
    [--java-opts "-Xms1g -Xmx8g"] \
    [--out-dir DIR] \
    [--primary-only] \
    --sample sample1.vcf.gz [--sample sample2.vcf.gz ...]

Purpose:
  Liftover one or more hs1/CHM13 VCF/gVCF files to GRCh38 for downstream tools
  (e.g. Exomiser/PharmCAT-compatible GRCh38 workflows).
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

CHAIN=""
REF_HG38=""
REF_ROOT="${REF_ROOT:-}"
PICARD_JAR=""
OUT_DIR=""
PRIMARY_ONLY=0
SAMPLES=""
JAVA_OPTS=""

while [ $# -gt 0 ]; do
  case "$1" in
    --chain-hs1-to-hg38) CHAIN="$2"; shift 2 ;;
    --ref-hg38) REF_HG38="$2"; shift 2 ;;
    --ref-root) REF_ROOT="$2"; shift 2 ;;
    --picard-jar) PICARD_JAR="$2"; shift 2 ;;
    --java-opts) JAVA_OPTS="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --primary-only) PRIMARY_ONLY=1; shift ;;
    --sample) SAMPLES="${SAMPLES} $2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$CHAIN" ] || { echo "ERROR: --chain-hs1-to-hg38 is required" >&2; exit 1; }
[ -n "$SAMPLES" ] || { echo "ERROR: at least one --sample is required" >&2; exit 1; }
[ -f "$CHAIN" ] || { echo "ERROR: not found: $CHAIN" >&2; exit 1; }

need_cmd java
need_cmd bcftools
need_cmd tabix
need_cmd samtools

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "${SCRIPT_DIR}/lib_liftover.sh"

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

if [ -z "$REF_ROOT" ]; then
  REF_ROOT="${BIO_TOOLS_INDEX_DIR:-${INDEX_DIR:-}}"
fi
if [ -z "$REF_HG38" ]; then
  roots=""
  [ -n "$REF_ROOT" ] && roots="$roots $REF_ROOT"
  roots="$roots ${PWD}/resources/ref ${SCRIPT_DIR}/../mapping/index"
  REF_HG38="$(pick_ref "$roots" \
    "GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna" \
    "hg38p14.H3630.fa.gz" \
    "hg38p14.H3630.fa" \
    "hg38p14DH3630O1102.fa" || true)"
fi
[ -n "$REF_HG38" ] || { echo "ERROR: --ref-hg38 not set and auto-detect failed." >&2; exit 1; }
[ -f "$REF_HG38" ] || { echo "ERROR: not found: $REF_HG38" >&2; exit 1; }

if [ -z "$PICARD_JAR" ]; then
  PICARD_JAR="$(detect_picard_jar "$SCRIPT_DIR" || true)"
fi
[ -n "$PICARD_JAR" ] || { echo "ERROR: set --picard-jar (picard.jar not auto-detected)" >&2; exit 1; }
[ -f "$PICARD_JAR" ] || { echo "ERROR: not found: $PICARD_JAR" >&2; exit 1; }
if [ -z "$JAVA_OPTS" ]; then
  JAVA_OPTS="${PICARD_JAVA_OPTS:-"-Xms1g -Xmx8g"}"
fi

ensure_ref_sidecars "$REF_HG38" "$PICARD_JAR" "$JAVA_OPTS"

if [ -z "$OUT_DIR" ]; then
  OUT_DIR="$(pwd)/experiments/liftover_hg38_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUT_DIR"

PRIMARY_CONTIGS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

REPORT="${OUT_DIR}/liftover_report.tsv"
echo "sample\tinput_type\tlifted_vcf\treject_vcf\tlifted_records\trejected_records\tprimary_records" > "$REPORT"

for sample in $SAMPLES; do
  [ -f "$sample" ] || { echo "ERROR: not found: $sample" >&2; exit 1; }

  in_type="vcf"
  case "$sample" in
    *.gvcf.gz|*.gvcf|*.g.vcf.gz|*.g.vcf) in_type="gvcf" ;;
  esac

  b=$(basename "$sample")
  b=${b%.g.vcf.gz}
  b=${b%.vcf.gz}
  b=${b%.gvcf.gz}
  b=${b%.g.vcf}
  b=${b%.vcf}
  b=${b%.gvcf}

  if [ "$in_type" = "gvcf" ]; then
    out_vcf="${OUT_DIR}/${b}.hg38.gvcf.gz"
    rej_vcf="${OUT_DIR}/${b}.hg38.gvcf.reject.vcf.gz"
  else
    out_vcf="${OUT_DIR}/${b}.hg38.vcf.gz"
    rej_vcf="${OUT_DIR}/${b}.hg38.reject.vcf.gz"
  fi

  [ ! -e "$out_vcf" ] || { echo "ERROR: output exists (would overwrite): $out_vcf" >&2; exit 2; }
  [ ! -e "$rej_vcf" ] || { echo "ERROR: output exists (would overwrite): $rej_vcf" >&2; exit 2; }

  echo "Liftover: $sample -> $out_vcf"
  java $JAVA_OPTS -jar "$PICARD_JAR" LiftoverVcf \
    I="$sample" \
    O="$out_vcf" \
    CHAIN="$CHAIN" \
    R="$REF_HG38" \
    REJECT="$rej_vcf" \
    RECOVER_SWAPPED_REF_ALT=true \
    WARN_ON_MISSING_CONTIG=true
  tabix -f "$out_vcf"
  tabix -f "$rej_vcf"

  if [ "$PRIMARY_ONLY" -eq 1 ]; then
    if [ "$in_type" = "gvcf" ]; then
      primary_vcf="${OUT_DIR}/${b}.hg38.gvcf.primary.vcf.gz"
    else
      primary_vcf="${OUT_DIR}/${b}.hg38.primary.vcf.gz"
    fi
    [ ! -e "$primary_vcf" ] || { echo "ERROR: output exists (would overwrite): $primary_vcf" >&2; exit 2; }
    bcftools view -r "$PRIMARY_CONTIGS" -Oz -o "$primary_vcf" "$out_vcf"
    tabix -f "$primary_vcf"
    primary_count=$(bcftools view -H "$primary_vcf" | wc -l)
  else
    primary_count="NA"
  fi

  lifted_count=$(bcftools view -H "$out_vcf" | wc -l)
  reject_count=$(bcftools view -H "$rej_vcf" | wc -l)

  echo "${sample}\t${in_type}\t${out_vcf}\t${rej_vcf}\t${lifted_count}\t${reject_count}\t${primary_count}" >> "$REPORT"
done

echo "Done. Report: $REPORT"
