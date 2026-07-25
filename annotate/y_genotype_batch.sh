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

set -uo pipefail

BAMDIR="${1:?usage: y_genotype_batch.sh <bam_dir> <out_dir> [marker_dir]}"
OUTDIR="${2:?out_dir required}"
MARKERDIR="${3:-markers}"

REF="${REF:?set REF to the reference the BAMs were mapped to}"
SITES="${SITES:-${MARKERDIR}/iceman_novel_candidates_usable8.tsv}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

shopt -s nullglob
BAMS=("${BAMDIR}"/*.rmdup.bam)
[ "${#BAMS[@]}" -gt 0 ] || { echo "no *.rmdup.bam in ${BAMDIR}" >&2; exit 1; }
SETS=("${MARKERDIR}"/*.txt)
[ "${#SETS[@]}" -gt 0 ] || { echo "no marker sets in ${MARKERDIR}" >&2; exit 1; }

mkdir -p "$OUTDIR"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "=== genotyping ${#BAMS[@]} sample(s) against ${#SETS[@]} marker set(s) ==="
echo "ref:   ${REF}"
echo "sites: ${SITES}"

# --- named marker sets -------------------------------------------------------
for set_file in "${SETS[@]}"; do
  set_name="$(basename "$set_file" .txt)"
  out="${OUTDIR}/swiss_${set_name}.tsv"
  : > "${TMP}/body"
  hdr=""

  for bam in "${BAMS[@]}"; do
    sample="$(basename "$bam" .rmdup.bam)"
    if ! "${HERE}/y_markers_pileup.py" --bam "$bam" --ref "$REF" \
           --marker-file "$set_file" --label "$set_name" \
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
  out="${OUTDIR}/swiss_novel_sites.tsv"
  : > "${TMP}/body"
  hdr=""
  for bam in "${BAMS[@]}"; do
    sample="$(basename "$bam" .rmdup.bam)"
    if ! "${HERE}/y_sites_pileup.py" --bam "$bam" --ref "$REF" \
           --sites "$SITES" --sample "$sample" > "${TMP}/one.tsv" 2>"${TMP}/err"; then
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
COV="${OUTDIR}/swiss_coverage.tsv"
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
echo "=== done ==="
