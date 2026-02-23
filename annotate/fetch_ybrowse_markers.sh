#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  fetch_ybrowse_markers.sh [--out PATH] [--refresh] [--fallback PATH]

Fetches YBrowse hg38 chrY marker VCF and index.
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
REFRESH=0
FALLBACK="${YMARKERS_FALLBACK:-}"
URL_BASE="http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz"

while [ $# -gt 0 ]; do
  case "$1" in
    --out) OUT="$2"; shift 2 ;;
    --refresh) REFRESH=1; shift ;;
    --fallback) FALLBACK="$2"; shift 2 ;;
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

echo "Marker header:"
bcftools view -h "$OUT" | awk '/^##date=|^##copyright=|^##contig=<ID=chrY/'
echo "Marker file ready: $OUT"
