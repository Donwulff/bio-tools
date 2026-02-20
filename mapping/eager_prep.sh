#!/bin/sh
set -eu

# Example helper for EAGER-style prefix repair + pairing.
# Override defaults via env vars for ad-hoc testing:
#   UBAM=/path/in.bam OUT=/path/out.bam THREADS=8 TMP_PREFIX=tmpfix ./eager_prep.sh

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)
EAGER_REPAIR="${REPO_ROOT}/util/eager_repair_bam.py"

UBAM="${UBAM:-/mnt/AncientDNA/Iceman-2024/iceman.oetzi.UDG_D2050_combined.mapped_rmdup.bam}"
OUT="${OUT:-iceman.oetzi.UDG_D2050_combined.mapped_rmdup.pair.bam}"
THREADS="${THREADS:-7}"
SORT_MEM="${SORT_MEM:-2G}"
TMP_PREFIX="${TMP_PREFIX:-foob49}"

samtools view -u -h "${UBAM}" \
  | python3 "${EAGER_REPAIR}" \
      - --strip-prefix -o - --output-uncompressed \
  | samtools sort -n -T "${TMP_PREFIX}" -m"${SORT_MEM}" -@"${THREADS}" -u \
  | python3 "${EAGER_REPAIR}" \
      - -o - --output-uncompressed \
      --pair --duplicate-merged --dup-orientation flip \
  | samtools fixmate -m - "${OUT}"
