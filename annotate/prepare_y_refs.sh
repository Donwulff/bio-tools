#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  prepare_y_refs.sh \
    --ref-hs1-src /path/to/hs1.fa \
    --ref-hg38-src /path/to/hg38.fa \
    [--out-dir ./resources/ref] \
    [--copy] \
    [--primary-only] \
    [--primary-contigs "chr1 ... chr22 chrX chrY chrM"] \
    [--hs1-name chm13v2.0_maskedY_rCRS.fa] \
    [--hg38-name hg38p14.H3630.fa]

Purpose:
  Stage hs1 and hg38 reference FASTAs for Y-analysis workflows and ensure:
  - FASTA index (.fai)
  - sequence dictionary (.dict)

Defaults:
  --out-dir defaults to ./resources/ref
  references are symlinked by default (use --copy to copy files)
  dictionary creation prefers: gatk, auto-detected picard.jar, picard CLI, or PICARD_JAR env
  --primary-only extracts canonical contigs from source FASTA into staged FASTA
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

REF_HS1_SRC=""
REF_HG38_SRC=""
OUT_DIR="${YREF_OUT_DIR:-$(pwd)/resources/ref}"
HS1_NAME="${YREF_HS1_NAME:-chm13v2.0_maskedY_rCRS.fa}"
HG38_NAME="${YREF_HG38_NAME:-hg38p14.H3630.fa}"
DO_COPY=0
PRIMARY_ONLY=0
PRIMARY_CONTIGS="${YREF_PRIMARY_CONTIGS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM}"
PICARD_JAR="${PICARD_JAR:-}"

while [ $# -gt 0 ]; do
  case "$1" in
    --ref-hs1-src) REF_HS1_SRC="$2"; shift 2 ;;
    --ref-hg38-src) REF_HG38_SRC="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --hs1-name) HS1_NAME="$2"; shift 2 ;;
    --hg38-name) HG38_NAME="$2"; shift 2 ;;
    --copy) DO_COPY=1; shift ;;
    --primary-only) PRIMARY_ONLY=1; shift ;;
    --primary-contigs) PRIMARY_CONTIGS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$REF_HS1_SRC" ] || { echo "ERROR: --ref-hs1-src is required" >&2; exit 1; }
[ -n "$REF_HG38_SRC" ] || { echo "ERROR: --ref-hg38-src is required" >&2; exit 1; }
[ -f "$REF_HS1_SRC" ] || { echo "ERROR: source not found: $REF_HS1_SRC" >&2; exit 1; }
[ -f "$REF_HG38_SRC" ] || { echo "ERROR: source not found: $REF_HG38_SRC" >&2; exit 1; }

need_cmd samtools
need_cmd java

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "${SCRIPT_DIR}/lib_liftover.sh"
if [ -z "$PICARD_JAR" ]; then
  PICARD_JAR="$(detect_picard_jar "$SCRIPT_DIR" || true)"
fi

mkdir -p "$OUT_DIR"

stage_ref() {
  stage_src="$1"
  stage_dst="$2"
  if [ "$PRIMARY_ONLY" -eq 1 ]; then
    # Build a staged FASTA containing only canonical contigs present in source.
    [ -f "${stage_src}.fai" ] || samtools faidx "$stage_src"
    keep=""
    for c in $PRIMARY_CONTIGS; do
      if awk -v cc="$c" '$1==cc{found=1} END{exit found?0:1}' "${stage_src}.fai"; then
        keep="$keep $c"
      fi
    done
    [ -n "$keep" ] || {
      echo "ERROR: none of the requested primary contigs were found in $stage_src" >&2
      exit 2
    }
    # shellcheck disable=SC2086
    samtools faidx "$stage_src" $keep > "$stage_dst"
    return 0
  fi
  if [ "$DO_COPY" -eq 1 ]; then
    cp -f "$stage_src" "$stage_dst"
  else
    ln -sfn "$stage_src" "$stage_dst"
  fi
}

copy_or_link_sidecar() {
  side_src="$1"
  side_dst="$2"
  if [ ! -f "$side_src" ]; then
    return 0
  fi
  if [ "$DO_COPY" -eq 1 ]; then
    cp -f "$side_src" "$side_dst"
  else
    ln -sfn "$side_src" "$side_dst"
  fi
}

ensure_dict() {
  ref="$1"
  ensure_ref_sidecars "$ref" "$PICARD_JAR" ""
}

prepare_one() {
  ref_src="$1"
  ref_dst="$2"
  stage_ref "$ref_src" "$ref_dst"

  if [ "$PRIMARY_ONLY" -eq 0 ]; then
    # Try to reuse existing sidecars from source location.
    copy_or_link_sidecar "${ref_src}.fai" "${ref_dst}.fai"
    copy_or_link_sidecar "${ref_src%.*}.dict" "${ref_dst%.*}.dict"
    copy_or_link_sidecar "${ref_src}.dict" "${ref_dst%.*}.dict"
  fi

  ensure_dict "$ref_dst"
}

HS1_DST="${OUT_DIR}/${HS1_NAME}"
HG38_DST="${OUT_DIR}/${HG38_NAME}"

prepare_one "$REF_HS1_SRC" "$HS1_DST"
prepare_one "$REF_HG38_SRC" "$HG38_DST"

echo "Prepared references:"
echo "  hs1  = $HS1_DST"
echo "  hg38 = $HG38_DST"
echo "  hs1 fai  = ${HS1_DST}.fai"
echo "  hg38 fai = ${HG38_DST}.fai"
echo "  hs1 dict  = $(dict_path_for_ref "$HS1_DST")"
echo "  hg38 dict = $(dict_path_for_ref "$HG38_DST")"
