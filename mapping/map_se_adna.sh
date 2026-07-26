#!/bin/bash
# Map single-end ancient DNA reads (short, variable length) to a GRCh38 build.
#
# The de-EAGER path in METHODS.md exists to undo read merging in already-mapped,
# pre-filtered BAMs. It does not apply to raw single-end FASTQ from SRA/ENA:
# there are no pairs to reconstruct and no prefixes to restore. This script is
# the entry point for that case.
#
# Aligner choice: bwa aln/samse, not bwa mem. bwa mem is designed for reads
# >=70 bp and loses sensitivity below that; typical aDNA capture reads are
# 35-90 bp. Defaults follow nf-core/eager 2.5.0 (-n 0.01 -k 2 -l 1024).
# -l exceeds any realistic aDNA read length, which disables seeding so the whole
# read -- including damage-bearing termini -- participates in the alignment.
#
# Caveat recorded deliberately: bwa aln is NOT alt-aware. bwa-postalt.js is not
# applied and the reference's .alt companion goes unused. Acceptable for Y
# marker genotyping, where MAPQ filtering does the equivalent work, but it is a
# real difference from any bwa mem run against the same reference.
#
# Deduplication is mandatory, not optional. Capture data is PCR-amplified; any
# downstream rule that counts "independent reads" is meaningless if duplicate
# copies of one molecule are allowed to inflate apparent support.
#
# A sample split across several sequencing runs is mapped run-by-run and merged
# BEFORE deduplication, never after. PCR duplicates of one library are duplicates
# whichever run they were sequenced in, so per-run markdup would leave
# cross-run copies of the same molecule in place and every "independent reads"
# threshold downstream would be inflated by exactly the reads it is meant to
# exclude. Each run keeps its own @RG, so the merge stays auditable.
#
# Usage:
#   mapping/map_se_adna.sh <fastq.gz> <sample_id> [run_id] [threads] [more.fastq.gz[:run_id] ...]
#
# The first four arguments are unchanged, so invocations recorded in RUNLOG.md
# before multi-run support existed still run identically. Additional runs for the
# same sample are appended as trailing arguments, optionally suffixed ":run_id"
# (default: the filename's <run> field, per the <alias>.<run>.fastq.gz convention).
#
# Environment:
#   REF      reference FASTA with a bwa index (required)
#   OUTDIR   output directory (default: current directory)
#   LOGDIR   log directory   (default: $OUTDIR/logs)
#   ALN_OPTS override bwa aln parameters
#   CALLABLE chrY callable/MSY denominator for the DoC figure (default 12570000)

set -euo pipefail

FQ="${1:?usage: map_se_adna.sh <fastq.gz> <sample_id> [run_id] [threads] [more.fastq.gz[:run] ...]}"
SAMPLE="${2:?sample_id required}"
RUN="${3:-$SAMPLE}"
THREADS="${4:-4}"
shift 4 2>/dev/null || shift $#

REF="${REF:?set REF to an indexed reference FASTA}"
OUTDIR="${OUTDIR:-.}"
LOGDIR="${LOGDIR:-$OUTDIR/logs}"
ALN_OPTS="${ALN_OPTS:--n 0.01 -k 2 -l 1024}"
CALLABLE="${CALLABLE:-12570000}"

for tool in bwa samtools; do
  command -v "$tool" >/dev/null || { echo "missing dependency: $tool" >&2; exit 1; }
done
[ -e "${REF}.bwt" ] || { echo "no bwa index beside $REF" >&2; exit 1; }

mkdir -p "$OUTDIR" "$LOGDIR"

# (fastq, run) work list: the positional pair first, then any trailing runs.
FQS=("$FQ");  RUNS=("$RUN")
for spec in "$@"; do
  case "$spec" in
    *.fastq.gz:*|*.fq.gz:*) p="${spec%:*}"; r="${spec##*:}" ;;
    *)                      p="$spec"
                            b="$(basename "$p")"; b="${b%.fastq.gz}"; b="${b%.fq.gz}"
                            r="${b#*.}"; [ "$r" = "$b" ] && r="$b" ;;
  esac
  [ -s "$p" ] || { echo "no such fastq: $p" >&2; exit 1; }
  FQS+=("$p"); RUNS+=("$r")
done

echo "=== ${SAMPLE} (${#FQS[@]} run(s): ${RUNS[*]}) threads=${THREADS} ==="
echo "ref:  ${REF}"
echo "opts: bwa aln ${ALN_OPTS}"

TIME=/usr/bin/time
[ -x "$TIME" ] || TIME=""

PARTS=()
for i in "${!FQS[@]}"; do
  fq="${FQS[$i]}"; run="${RUNS[$i]}"
  tag="${SAMPLE}"; [ "${#FQS[@]}" -gt 1 ] && tag="${SAMPLE}.${run}"
  RG="@RG\tID:${run}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPU:${run}\tPL:ILLUMINA"

  ${TIME:+$TIME -v} bwa aln -t "$THREADS" ${ALN_OPTS} "$REF" "$fq" \
      > "${OUTDIR}/${tag}.sai" 2> "${LOGDIR}/${tag}.aln.log"

  bwa samse -r "$RG" "$REF" "${OUTDIR}/${tag}.sai" "$fq" 2> "${LOGDIR}/${tag}.samse.log" \
    | samtools sort -@2 -m 1G -T "${OUTDIR}/tmp.${tag}" -o "${OUTDIR}/${tag}.sorted.bam" -
  rm -f "${OUTDIR}/${tag}.sai"
  PARTS+=("${OUTDIR}/${tag}.sorted.bam")
done

if [ "${#PARTS[@]}" -gt 1 ]; then
  # Merge before markdup -- see the header note on cross-run duplicates.
  samtools merge -f -@2 "${OUTDIR}/${SAMPLE}.sorted.bam" "${PARTS[@]}"
  rm -f "${PARTS[@]}"
fi
samtools index "${OUTDIR}/${SAMPLE}.sorted.bam"

samtools markdup -r -s "${OUTDIR}/${SAMPLE}.sorted.bam" "${OUTDIR}/${SAMPLE}.rmdup.bam" \
    2> "${LOGDIR}/${SAMPLE}.markdup.log"
samtools index "${OUTDIR}/${SAMPLE}.rmdup.bam"

BAM="${OUTDIR}/${SAMPLE}.rmdup.bam"
{
  echo "sample=${SAMPLE} runs=${RUNS[*]} n_runs=${#FQS[@]}"
  for i in "${!FQS[@]}"; do
    tag="${SAMPLE}"; [ "${#FQS[@]}" -gt 1 ] && tag="${SAMPLE}.${RUNS[$i]}"
    awk -v r="${RUNS[$i]}" '/Maximum resident/{printf "peak_rss_gb[%s]=%.2f\n",r,$NF/1048576}' \
        "${LOGDIR}/${tag}.aln.log" 2>/dev/null
    awk -v r="${RUNS[$i]}" '/wall clock/{print "aln_wall["r"]="$NF}' \
        "${LOGDIR}/${tag}.aln.log" 2>/dev/null
  done
  samtools idxstats "$BAM" | awk -v OFS== '$1=="chrY"{print "chrY_reads",$3}'
  echo "chrY_reads_mq25=$(samtools view -c -q 25 "$BAM" chrY)"
  samtools depth -a -q 20 -Q 25 -r chrY "$BAM" 2>/dev/null \
    | awk -v d="$CALLABLE" '{s+=$3} END{printf "chrY_bases_mq25_bq20=%d\nchrY_DoC_callable=%.4fx\n", s, s/d}'
} | tee "${LOGDIR}/${SAMPLE}.summary"

echo "=== ${SAMPLE} done ==="
