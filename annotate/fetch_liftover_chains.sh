#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  fetch_liftover_chains.sh [--out-dir DIR] [--refresh]

Fetches GRCh38<->hs1(CHM13) chain files using multiple mirrors.

Outputs (if successful):
  <out-dir>/hg38ToHs1.over.chain.gz
  <out-dir>/hs1ToHg38.over.chain.gz

Notes:
  - If UCSC hosts are unreachable, the script also tries MARBL/CHM13 chain names.
  - Existing files are reused unless --refresh is provided.
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

OUT_DIR="$(pwd)/resources/chains"
REFRESH=0

while [ $# -gt 0 ]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --refresh) REFRESH=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

need_cmd wget

mkdir -p "$OUT_DIR"

HG38_TO_HS1="${OUT_DIR}/hg38ToHs1.over.chain.gz"
HS1_TO_HG38="${OUT_DIR}/hs1ToHg38.over.chain.gz"

fetch_one() {
  out="$1"
  shift
  if [ -f "$out" ] && [ "$REFRESH" -eq 0 ]; then
    echo "Using existing: $out"
    return 0
  fi
  tmp="${out}.tmp"
  rm -f "$tmp"
  for u in "$@"; do
    echo "Try: $u"
    if wget --inet4-only -nv -O "$tmp" "$u"; then
      mv -f "$tmp" "$out"
      echo "Fetched: $out"
      return 0
    fi
  done
  rm -f "$tmp"
  return 1
}

ok1=0
ok2=0

if fetch_one "$HG38_TO_HS1" \
  "https://hgdownload.gi.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz" \
  "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain"; then
  ok1=1
fi

if fetch_one "$HS1_TO_HG38" \
  "https://hgdownload.gi.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz" \
  "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain"; then
  ok2=1
fi

if [ "$ok1" -ne 1 ] || [ "$ok2" -ne 1 ]; then
  echo "ERROR: could not fetch all required chain files." >&2
  echo "Expected outputs:" >&2
  echo "  $HG38_TO_HS1" >&2
  echo "  $HS1_TO_HG38" >&2
  exit 2
fi

echo "Chain files ready:"
echo "  $HG38_TO_HS1"
echo "  $HS1_TO_HG38"
