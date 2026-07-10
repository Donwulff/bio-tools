#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  fetch_ybrowse_markers.sh [--out PATH] [--gff3-out PATH] [--refresh] [--fallback PATH] [--no-gff3]

Fetches YBrowse hg38 chrY marker VCF (+index) and optionally GFF3.
Default output:
  ./resources/snps_hg38.vcf.gz
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

OUT="${YMARKERS_OUT:-$(pwd)/resources/snps_hg38.vcf.gz}"
GFF3_OUT="${YMARKERS_GFF3_OUT:-}"
REFRESH=0
FALLBACK="${YMARKERS_FALLBACK:-}"
FETCH_GFF3=1
URL_BASE="http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz"
URL_GFF3="http://ybrowse.org/gbrowse2/gff/snps_hg38.gff3"

# Optional local config (same convention as mapping scripts).
if [ -f "${BIO_TOOLS_CFG:-./bio-tools.cfg}" ]; then
  # shellcheck disable=SC1090
  . "${BIO_TOOLS_CFG:-./bio-tools.cfg}"
  if [ -z "${YMARKERS_OUT:-}" ] && [ -n "${Y_RESOURCES_DIR:-}" ]; then
    OUT="${Y_RESOURCES_DIR%/}/snps_hg38.vcf.gz"
  fi
fi

while [ $# -gt 0 ]; do
  case "$1" in
    --out) OUT="$2"; shift 2 ;;
    --gff3-out) GFF3_OUT="$2"; shift 2 ;;
    --refresh) REFRESH=1; shift ;;
    --fallback) FALLBACK="$2"; shift 2 ;;
    --no-gff3) FETCH_GFF3=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

need_cmd bcftools
need_cmd tabix
need_cmd wget
need_cmd sh

OUT_DIR=$(dirname "$OUT")
mkdir -p "$OUT_DIR"
if [ -z "$GFF3_OUT" ]; then
  case "$OUT" in
    *.vcf.gz) GFF3_OUT="${OUT%.vcf.gz}.gff3" ;;
    *) GFF3_OUT="${OUT}.gff3" ;;
  esac
fi

if [ -f "$OUT" ] && [ -f "${OUT}.tbi" ] && [ "$REFRESH" -eq 0 ]; then
  echo "Using existing marker file: $OUT"
else
  echo "Fetching latest markers from YBrowse..."
  set +e
  wget -nv -O "$OUT" "$URL_BASE"
  rc1=$?
  wget -nv -O "${OUT}.tbi" "${URL_BASE}.tbi"
  rc2=$?
  set -e
  if [ "$rc1" -ne 0 ] || [ "$rc2" -ne 0 ]; then
    echo "WARNING: download failed (rc=$rc1/$rc2)." >&2
    if [ -n "$FALLBACK" ] && [ -f "$FALLBACK" ] && [ -f "${FALLBACK}.tbi" ]; then
      echo "Using fallback marker file: $FALLBACK" >&2
      cp -f "$FALLBACK" "$OUT"
      cp -f "${FALLBACK}.tbi" "${OUT}.tbi"
    else
      echo "ERROR: no usable fallback markers available." >&2
      exit 2
    fi
  fi
fi

if [ ! -f "${OUT}.tbi" ]; then
  tabix -f "$OUT"
fi

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SAN_OUT="${OUT%.vcf.gz}.sanitized.vcf.gz"
"${SCRIPT_DIR}/sanitize_marker_vcf.sh" --input "$OUT" --output "$SAN_OUT"
mv -f "$SAN_OUT" "$OUT"
mv -f "${SAN_OUT}.tbi" "${OUT}.tbi"

if [ "$FETCH_GFF3" -eq 1 ]; then
  mkdir -p "$(dirname "$GFF3_OUT")"
  if [ -f "$GFF3_OUT" ] && [ "$REFRESH" -eq 0 ]; then
    echo "Using existing marker GFF3: $GFF3_OUT"
  else
    echo "Fetching latest marker GFF3 from YBrowse..."
    set +e
    wget -nv -O "$GFF3_OUT" "$URL_GFF3"
    rc3=$?
    set -e
    if [ "$rc3" -ne 0 ]; then
      echo "WARNING: GFF3 download failed (rc=$rc3)." >&2
      # If VCF fallback was used, try sibling .gff3 as convenience.
      if [ -n "$FALLBACK" ]; then
        fb_gff3="${FALLBACK%.vcf.gz}.gff3"
        if [ -f "$fb_gff3" ]; then
          echo "Using fallback marker GFF3: $fb_gff3" >&2
          cp -f "$fb_gff3" "$GFF3_OUT"
        fi
      fi
    fi
  fi
fi

echo "Marker header:"
bcftools view -h "$OUT" | awk '/^##date=|^##copyright=|^##contig=<ID=chrY/'
echo "Marker file ready: $OUT"
if [ "$FETCH_GFF3" -eq 1 ] && [ -f "$GFF3_OUT" ]; then
  echo "Marker GFF3 ready: $GFF3_OUT"
fi
