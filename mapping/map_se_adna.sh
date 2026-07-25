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
# Usage:
#   mapping/map_se_adna.sh <fastq.gz> <sample_id> [run_id] [threads]
# Environment:
#   REF      reference FASTA with a bwa index (required)
#   OUTDIR   output directory (default: current directory)
#   LOGDIR   log directory   (default: $OUTDIR/logs)
#   ALN_OPTS override bwa aln parameters
#   CALLABLE chrY callable/MSY denominator for the DoC figure (default 12570000)

set -euo pipefail

FQ="${1:?usage: map_se_adna.sh <fastq.gz> <sample_id> [run_id] [threads]}"
SAMPLE="${2:?sample_id required}"
RUN="${3:-$SAMPLE}"
THREADS="${4:-4}"

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
RG="@RG\tID:${RUN}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPU:${RUN}\tPL:ILLUMINA"

echo "=== ${SAMPLE} (${RUN}) threads=${THREADS} ==="
echo "ref:  ${REF}"
echo "opts: bwa aln ${ALN_OPTS}"

TIME=/usr/bin/time
[ -x "$TIME" ] || TIME=""

${TIME:+$TIME -v} bwa aln -t "$THREADS" ${ALN_OPTS} "$REF" "$FQ" \
    > "${OUTDIR}/${SAMPLE}.sai" 2> "${LOGDIR}/${SAMPLE}.aln.log"

bwa samse -r "$RG" "$REF" "${OUTDIR}/${SAMPLE}.sai" "$FQ" 2> "${LOGDIR}/${SAMPLE}.samse.log" \
  | samtools sort -@2 -m 1G -T "${OUTDIR}/tmp.${SAMPLE}" -o "${OUTDIR}/${SAMPLE}.sorted.bam" -
samtools index "${OUTDIR}/${SAMPLE}.sorted.bam"
rm -f "${OUTDIR}/${SAMPLE}.sai"

samtools markdup -r -s "${OUTDIR}/${SAMPLE}.sorted.bam" "${OUTDIR}/${SAMPLE}.rmdup.bam" \
    2> "${LOGDIR}/${SAMPLE}.markdup.log"
samtools index "${OUTDIR}/${SAMPLE}.rmdup.bam"

BAM="${OUTDIR}/${SAMPLE}.rmdup.bam"
{
  echo "sample=${SAMPLE} run=${RUN}"
  awk '/Maximum resident/{printf "peak_rss_gb=%.2f\n",$NF/1048576}' "${LOGDIR}/${SAMPLE}.aln.log" 2>/dev/null
  awk '/wall clock/{print "aln_wall="$NF}' "${LOGDIR}/${SAMPLE}.aln.log" 2>/dev/null
  samtools idxstats "$BAM" | awk -v OFS== '$1=="chrY"{print "chrY_reads",$3}'
  echo "chrY_reads_mq25=$(samtools view -c -q 25 "$BAM" chrY)"
  samtools depth -a -q 20 -Q 25 -r chrY "$BAM" 2>/dev/null \
    | awk -v d="$CALLABLE" '{s+=$3} END{printf "chrY_bases_mq25_bq20=%d\nchrY_DoC_callable=%.4fx\n", s, s/d}'
} | tee "${LOGDIR}/${SAMPLE}.summary"

echo "=== ${SAMPLE} done ==="
