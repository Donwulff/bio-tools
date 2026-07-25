#!/bin/bash
# Map a staging directory of single-end aDNA FASTQs, one sample at a time.
#
# Thin driver over mapping/map_se_adna.sh. It exists so that a mapping run is
# reproducible from the repository rather than from a shell loop typed once and
# lost -- the parameters below are the record of how the BAMs were made.
#
# Samples are processed SEQUENTIALLY on purpose. bwa aln holds the whole BWT+SA
# index resident (~4 GB against the hg38p14DH3630O build, measured peak RSS
# 4.05 GB at -t 4); running several at once multiplies that by the job count and
# has already caused an OOM kill on this machine. Parallelism belongs in -t,
# which shares one copy of the index.
#
# Resumable: a sample whose .rmdup.bam already exists with an index is skipped,
# so an interrupted run is restarted by re-invoking the same command.
#
# Input files are named <alias>.<run>.fastq.gz, the convention written by
# annotate/fetch_ena_runs.sh; alias becomes the sample ID and run the read group.
#
# Usage:
#   mapping/map_se_batch.sh <staging_dir> [threads]
# Environment:
#   REF      reference FASTA with a bwa index (required, passed through)
#   OUTDIR   output directory (default: <staging_dir>/bam)
#   ALN_OPTS override bwa aln parameters (passed through)

set -uo pipefail

STAGE="${1:?usage: map_se_batch.sh <staging_dir> [threads]}"
THREADS="${2:-8}"

REF="${REF:?set REF to an indexed reference FASTA}"
OUTDIR="${OUTDIR:-${STAGE}/bam}"
export REF OUTDIR
LOGDIR="${LOGDIR:-${OUTDIR}/logs}"
export LOGDIR

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAPPER="${HERE}/map_se_adna.sh"
[ -x "$MAPPER" ] || { echo "not executable: $MAPPER" >&2; exit 1; }

shopt -s nullglob
FQS=("${STAGE}"/*.fastq.gz)
[ "${#FQS[@]}" -gt 0 ] || { echo "no *.fastq.gz in ${STAGE}" >&2; exit 1; }

mkdir -p "$OUTDIR" "$LOGDIR"
echo "=== batch: ${#FQS[@]} fastq(s) from ${STAGE}, threads=${THREADS} ==="
echo "ref:    ${REF}"
echo "outdir: ${OUTDIR}"
echo "started $(date -Is)"

done_n=0 skip_n=0 fail_n=0
for fq in "${FQS[@]}"; do
  base="$(basename "$fq" .fastq.gz)"
  sample="${base%%.*}"                 # <alias> from <alias>.<run>
  run="${base#*.}"                     # <run>
  [ "$run" = "$base" ] && run="$sample"

  bam="${OUTDIR}/${sample}.rmdup.bam"
  if [ -s "$bam" ] && [ -s "${bam}.bai" ]; then
    echo "SKIP  ${sample} (${bam} present)"
    skip_n=$((skip_n + 1))
    continue
  fi

  echo
  echo "---- ${sample} [$((done_n + skip_n + fail_n + 1))/${#FQS[@]}] $(date -Is) ----"
  if "$MAPPER" "$fq" "$sample" "$run" "$THREADS"; then
    done_n=$((done_n + 1))
  else
    echo "FAIL  ${sample} (exit $?)" >&2
    rm -f "${OUTDIR}/${sample}.sai"    # truncated .sai would poison a rerun
    fail_n=$((fail_n + 1))
  fi
done

echo
echo "=== batch complete $(date -Is): mapped=${done_n} skipped=${skip_n} failed=${fail_n} ==="
[ "$fail_n" -eq 0 ]
