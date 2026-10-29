#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
. "${SCRIPT_DIR}/lib_liftover.sh"
load_liftover_config

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
    [--backend auto|gatk|picard-jar|picard-cli] \
    [--picard-jar picard.jar] [--java-opts "-Xms1g -Xmx8g"] \
    [--hs1-name chm13v2.0_maskedY_rCRS.fa] \
    [--hg38-name hg38.fa]

Purpose:
  Stage hs1 and hg38 reference FASTAs for Y-analysis workflows and ensure:
  - FASTA index (.fai)
  - sequence dictionary (.dict)

Defaults:
  --out-dir defaults to Y_RESOURCES_DIR/ref or ./resources/ref
  references are symlinked by default (use --copy to copy files)
  dictionary creation shares backend selection/JVM options with the liftover drivers
  --primary-only extracts canonical contigs from source FASTA into staged FASTA
  existing targets are refused except unchanged source symlinks; use a fresh out-dir
  validation requires python3 + pysam; ordinary gzip must be staged uncompressed or BGZF
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
OUT_DIR="${YREF_OUT_DIR:-${Y_RESOURCES_DIR:-$(pwd)/resources}/ref}"
HS1_NAME="${YREF_HS1_NAME:-chm13v2.0_maskedY_rCRS.fa}"
HG38_NAME="${YREF_HG38_NAME:-hg38.fa}"
DO_COPY=0
PRIMARY_ONLY=0
PRIMARY_CONTIGS="${YREF_PRIMARY_CONTIGS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM}"
PICARD_JAR="${PICARD_JAR:-}"
LIFTOVER_BACKEND="${LIFTOVER_BACKEND:-auto}"
LIFTOVER_JAVA_OPTS="${PICARD_JAVA_OPTS-"-Xms1g -Xmx8g"}"

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
    --backend) LIFTOVER_BACKEND="$2"; shift 2 ;;
    --picard-jar) PICARD_JAR="$2"; shift 2 ;;
    --java-opts) LIFTOVER_JAVA_OPTS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$REF_HS1_SRC" ] || { echo "ERROR: --ref-hs1-src is required" >&2; exit 1; }
[ -n "$REF_HG38_SRC" ] || { echo "ERROR: --ref-hg38-src is required" >&2; exit 1; }
[ -f "$REF_HS1_SRC" ] || { echo "ERROR: source not found: $REF_HS1_SRC" >&2; exit 1; }
[ -f "$REF_HG38_SRC" ] || { echo "ERROR: source not found: $REF_HG38_SRC" >&2; exit 1; }

need_cmd samtools
need_cmd python3

python3 "${SCRIPT_DIR}/liftover_backend.py" --check-options \
  --backend "$LIFTOVER_BACKEND" --java-opts="$LIFTOVER_JAVA_OPTS"

OUT_DIR=$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().absolute())' "$OUT_DIR")
REF_HS1_SRC=$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).absolute())' "$REF_HS1_SRC")
REF_HG38_SRC=$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).absolute())' "$REF_HG38_SRC")

check_stage_target() {
  if [ -e "$2" ] || [ -L "$2" ]; then
    if [ "$PRIMARY_ONLY" -eq 0 ] && [ "$DO_COPY" -eq 0 ] && [ -L "$2" ] && [ "$1" -ef "$2" ]; then
      return 0
    fi
    echo "ERROR: refusing to replace staged reference $2; use a fresh --out-dir." >&2
    exit 2
  fi
  for side in "$2.fai" "$2.gzi" "$(dict_path_for_ref "$2")" "${2%.*}.dict"; do
    if [ ! -e "$2" ] && [ ! -L "$2" ] && { [ -e "$side" ] || [ -L "$side" ]; }; then
      echo "ERROR: orphaned sidecar $side; use a fresh --out-dir." >&2
      exit 2
    fi
  done
}

stage_ref() {
  stage_src="$1"
  stage_dst="$2"
  ensure_ref_index "$stage_src"
  if [ "$PRIMARY_ONLY" -eq 1 ]; then
    # Build a staged FASTA containing only canonical contigs present in source.
    set --
    for c in $PRIMARY_CONTIGS; do
      if awk -v cc="$c" '$1==cc{found=1} END{exit found?0:1}' "${stage_src}.fai"; then
        set -- "$@" "$c"
      fi
    done
    [ "$#" -gt 0 ] || {
      echo "ERROR: none of the requested primary contigs were found in $stage_src" >&2
      exit 2
    }
    stage_tmp=$(mktemp "${stage_dst}.tmp.XXXXXX")
    if samtools faidx "$stage_src" "$@" > "$stage_tmp"; then
      mv "$stage_tmp" "$stage_dst"
    else
      rm -f "$stage_tmp"
      return 2
    fi
    return 0
  fi
  if [ "$DO_COPY" -eq 1 ]; then
    cp "$stage_src" "$stage_dst"
  else
    [ -L "$stage_dst" ] || ln -s "$stage_src" "$stage_dst"
  fi
}

copy_or_link_sidecar() {
  side_src="$1"
  side_dst="$2"
  if [ ! -f "$side_src" ]; then
    return 0
  fi
  if [ -e "$side_dst" ] || [ -L "$side_dst" ]; then
    return 0
  fi
  if [ "$DO_COPY" -eq 1 ]; then
    cp "$side_src" "$side_dst"
  else
    ln -s "$side_src" "$side_dst"
  fi
}

ensure_dict() {
  ref="$1"
  # Reusing valid dictionaries does not require an installed Java backend.
  ensure_ref_sidecars "$ref" "$PICARD_JAR" "$LIFTOVER_JAVA_OPTS" "$BACKEND_SPEC"
}

prepare_one() {
  ref_src="$1"
  ref_dst="$2"
  stage_ref "$ref_src" "$ref_dst"

  if [ "$PRIMARY_ONLY" -eq 0 ]; then
    # Try to reuse existing sidecars from source location.
    copy_or_link_sidecar "${ref_src}.fai" "${ref_dst}.fai"
    copy_or_link_sidecar "${ref_src}.gzi" "${ref_dst}.gzi"
    copy_or_link_sidecar "$(dict_path_for_ref "$ref_src")" "$(dict_path_for_ref "$ref_dst")"
    copy_or_link_sidecar "${ref_src%.*}.dict" "$(dict_path_for_ref "$ref_dst")"
  fi

  ensure_dict "$ref_dst"
}

stage_path() {
  python3 -c 'from pathlib import Path; import sys; p=Path(sys.argv[1]); print(p.parent.resolve()/p.name)' "$1"
}

HS1_DST=$(stage_path "${OUT_DIR}/${HS1_NAME}")
HG38_DST=$(stage_path "${OUT_DIR}/${HG38_NAME}")
[ "$HS1_DST" != "$HG38_DST" ] || { echo "ERROR: reference output names must differ" >&2; exit 2; }
check_stage_target "$REF_HS1_SRC" "$HS1_DST"
check_stage_target "$REF_HG38_SRC" "$HG38_DST"

has_dictionary() {
  for candidate in "$(dict_path_for_ref "$1")" "${1%.*}.dict"; do
    if [ -f "$candidate" ]; then return 0; fi
  done
  return 1
}

# Select/validate transport before copying or subsetting whenever a dictionary
# will need generation. Cached dictionaries still need no installed backend.
BACKEND_SPEC=""
for reference in "$REF_HS1_SRC" "$REF_HG38_SRC" "$HS1_DST" "$HG38_DST"; do
  for candidate in "$(dict_path_for_ref "$reference")" "${reference%.*}.dict"; do
    if { [ -e "$candidate" ] || [ -L "$candidate" ]; } && [ ! -f "$candidate" ]; then
      echo "ERROR: invalid dictionary path (not a readable file target): $candidate" >&2
      exit 2
    fi
  done
done
if [ "$PRIMARY_ONLY" -eq 1 ] ||
   { ! has_dictionary "$REF_HS1_SRC" && ! has_dictionary "$HS1_DST"; } ||
   { ! has_dictionary "$REF_HG38_SRC" && ! has_dictionary "$HG38_DST"; }; then
  BACKEND_SPEC="$(select_liftover_backend "$LIFTOVER_BACKEND" "$PICARD_JAR" "$LIFTOVER_JAVA_OPTS")"
fi

mkdir -p "$OUT_DIR"
prepare_one "$REF_HS1_SRC" "$HS1_DST"
prepare_one "$REF_HG38_SRC" "$HG38_DST"

echo "Prepared references:"
echo "  hs1  = $HS1_DST"
echo "  hg38 = $HG38_DST"
echo "  hs1 fai  = ${HS1_DST}.fai"
echo "  hg38 fai = ${HG38_DST}.fai"
echo "  hs1 dict  = $(dict_path_for_ref "$HS1_DST")"
echo "  hg38 dict = $(dict_path_for_ref "$HG38_DST")"
