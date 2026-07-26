#!/bin/bash
# Genotype a directory of BAMs against the committed marker and site sets.
#
# Runs the two pileup tools across every sample and collates the results into
# one table per marker set, with a sample column that neither tool emits on its
# own. This is the step that turns mapped BAMs into the evidence tables a
# write-up cites, so it is a script rather than a loop typed at a prompt.
#
# Named markers go through annotate/y_markers_pileup.py, which resolves them via
# resources/marker_index.tsv.gz. Uncatalogued positions cannot be resolved by
# name -- being uncatalogued is what makes them novel -- so they go through
# annotate/y_sites_pileup.py by coordinate instead.
#
# Nothing here assigns a haplogroup. Output is per-SNP derived/ancestral/nocall
# by design; collapsing that into a terminal label is the error class this
# repository exists to document.
#
# Usage:
#   annotate/y_genotype_batch.sh <bam_dir> <out_dir> [marker_dir]
# Environment:
#   REF      reference FASTA the BAMs were mapped to (required)
#   SITES    site TSV for the by-coordinate pass
#            (default: <marker_dir>/iceman_novel_candidates_usable8.tsv)
#   PREFIX   output filename prefix (default: y)
#   INDEX    marker coordinate index (default: resources/marker_index.tsv.gz)
#   MINMQ    minimum MAPQ  (default 25, passed to both pileup tools)
#   MINBQ    minimum BASEQ (default 20, passed to both pileup tools)
#   MAXMQ0   pct_mq0 rejection threshold (default 30.0, via Y_MAX_PCT_MQ0)
#
# The filter parameters are injectable because a filtering criterion that turns
# out to be wrong should be re-runnable in one command rather than requiring an
# edit. MAPQ in particular is known to be allele-biased at some loci (NOTES,
# 2026-07-26): a uniquely-mapped derived read can fall below 25 purely because a
# suboptimal hit exists elsewhere. Every run writes <prefix>_params.txt so any
# table in results/ can be traced to the exact thresholds that produced it.
#
# PREFIX exists because this prefix was once hardcoded to "swiss", which would
# label Hungarian, Sardinian and Bavarian results as Swiss. The tables already
# committed under results/swiss/ and results/swiss15/ were produced with that
# hardcoded value and are reproduced by passing PREFIX=swiss explicitly.
#
# Not every .txt in the marker directory is a marker set -- markers/ also holds
# membership lists such as family_a_members.txt, whose entries are sample IDs.
# Globbing those in would genotype "MX150" as though it were a marker name, so
# each candidate file is checked against the marker index first and skipped, out
# loud, if none of its entries resolve.

set -uo pipefail

BAMDIR="${1:?usage: y_genotype_batch.sh <bam_dir> <out_dir> [marker_dir]}"
OUTDIR="${2:?out_dir required}"
MARKERDIR="${3:-markers}"

REF="${REF:?set REF to the reference the BAMs were mapped to}"
SITES="${SITES:-${MARKERDIR}/iceman_novel_candidates_usable8.tsv}"
PREFIX="${PREFIX:-y}"
INDEX="${INDEX:-resources/marker_index.tsv.gz}"
MINMQ="${MINMQ:-25}"
MINBQ="${MINBQ:-20}"
MAXMQ0="${MAXMQ0:-30.0}"
export Y_MAX_PCT_MQ0="$MAXMQ0"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

shopt -s nullglob
BAMS=("${BAMDIR}"/*.rmdup.bam)
[ "${#BAMS[@]}" -gt 0 ] || { echo "no *.rmdup.bam in ${BAMDIR}" >&2; exit 1; }

CAND=("${MARKERDIR}"/*.txt)
[ "${#CAND[@]}" -gt 0 ] || { echo "no .txt files in ${MARKERDIR}" >&2; exit 1; }

[ -s "$INDEX" ] || { echo "no marker index at ${INDEX}" >&2; exit 1; }
KNOWN="$(mktemp)" || exit 1
# LC_ALL=C throughout: sort and comm must agree on collation or comm silently
# mis-reports matches, and this comparison decides whether a file is treated as
# a marker set at all. Under a non-C locale it warns "input is not in sorted
# order" and can undercount, which would skip a real marker set.
zcat -f "$INDEX" | tail -n +2 | cut -f1 | LC_ALL=C sort -u > "$KNOWN"

SETS=()
for f in "${CAND[@]}"; do
  n=$(sed 's/#.*//' "$f" | tr -d ' \t' | grep -v '^$' | LC_ALL=C sort -u \
        | LC_ALL=C comm -12 - "$KNOWN" | wc -l)
  if [ "$n" -gt 0 ]; then
    SETS+=("$f")
  else
    echo "SKIP  $(basename "$f"): no entries resolve as markers (not a marker set)" >&2
  fi
done
rm -f "$KNOWN"
[ "${#SETS[@]}" -gt 0 ] || { echo "no usable marker sets in ${MARKERDIR}" >&2; exit 1; }

mkdir -p "$OUTDIR"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "=== genotyping ${#BAMS[@]} sample(s) against ${#SETS[@]} marker set(s) ==="
echo "ref:   ${REF}"
echo "sites: ${SITES}"

# --- named marker sets -------------------------------------------------------
for set_file in "${SETS[@]}"; do
  set_name="$(basename "$set_file" .txt)"
  out="${OUTDIR}/${PREFIX}_${set_name}.tsv"
  : > "${TMP}/body"
  hdr=""

  for bam in "${BAMS[@]}"; do
    sample="$(basename "$bam" .rmdup.bam)"
    if ! "${HERE}/y_markers_pileup.py" --bam "$bam" --ref "$REF" \
           --marker-file "$set_file" --label "$set_name" \
           --min-mq "$MINMQ" --min-bq "$MINBQ" \
           --out "${TMP}/one.tsv" >/dev/null 2>"${TMP}/err"; then
      echo "FAIL  ${sample} / ${set_name}: $(tail -1 "${TMP}/err")" >&2
      continue
    fi
    [ -n "$hdr" ] || hdr="sample	$(head -1 "${TMP}/one.tsv")"
    tail -n +2 "${TMP}/one.tsv" | sed "s/^/${sample}\t/" >> "${TMP}/body"
  done

  { echo "$hdr"; cat "${TMP}/body"; } > "$out"
  echo "  ${out}  ($(( $(wc -l < "$out") - 1 )) rows)"
done

# --- uncatalogued positions --------------------------------------------------
if [ -s "$SITES" ]; then
  out="${OUTDIR}/${PREFIX}_novel_sites.tsv"
  : > "${TMP}/body"
  hdr=""
  for bam in "${BAMS[@]}"; do
    sample="$(basename "$bam" .rmdup.bam)"
    if ! "${HERE}/y_sites_pileup.py" --bam "$bam" --ref "$REF" \
           --sites "$SITES" --sample "$sample" \
           --min-mq "$MINMQ" --min-bq "$MINBQ" > "${TMP}/one.tsv" 2>"${TMP}/err"; then
      echo "FAIL  ${sample} / novel sites: $(tail -1 "${TMP}/err")" >&2
      continue
    fi
    [ -n "$hdr" ] || hdr="$(head -1 "${TMP}/one.tsv")"
    tail -n +2 "${TMP}/one.tsv" >> "${TMP}/body"
  done
  { echo "$hdr"; cat "${TMP}/body"; } > "$out"
  echo "  ${out}  ($(( $(wc -l < "$out") - 1 )) rows)"
fi

# --- coverage, always with its denominator -----------------------------------
COV="${OUTDIR}/${PREFIX}_coverage.tsv"
{
  printf 'sample\tchrY_reads\tchrY_reads_mq25\tchrY_bases_mq25_bq20\tchrY_DoC_callable\n'
  for bam in "${BAMS[@]}"; do
    sample="$(basename "$bam" .rmdup.bam)"
    s="${BAMDIR}/logs/${sample}.summary"
    [ -s "$s" ] || { printf '%s\tNA\tNA\tNA\tNA\n' "$sample"; continue; }
    awk -F= -v s="$sample" '
      $1=="chrY_reads"{a=$2} $1=="chrY_reads_mq25"{b=$2}
      $1=="chrY_bases_mq25_bq20"{c=$2} $1=="chrY_DoC_callable"{d=$2}
      END{printf "%s\t%s\t%s\t%s\t%s\n", s, a?a:"NA", b?b:"NA", c?c:"NA", d?d:"NA"}' "$s"
  done
} > "$COV"
echo "  ${COV}"

# Provenance: a results table is only reproducible if it records what produced
# it. Written last so it is absent if the run failed part-way.
PAR="${OUTDIR}/${PREFIX}_params.txt"
{
  echo "generated=$(date -Is)"
  echo "bam_dir=${BAMDIR}"
  echo "marker_dir=${MARKERDIR}"
  echo "ref=${REF}"
  echo "sites=${SITES}"
  echo "index=${INDEX}"
  echo "min_mq=${MINMQ}"
  echo "min_bq=${MINBQ}"
  echo "max_pct_mq0=${MAXMQ0}"
  echo "n_samples=${#BAMS[@]}"
  echo "marker_sets=$(printf '%s ' "${SETS[@]##*/}")"
  echo "git_commit=$(git -C "$HERE" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "git_dirty=$(git -C "$HERE" status --porcelain 2>/dev/null | wc -l)"
} > "$PAR"
echo "  ${PAR}"
echo "=== done ==="
