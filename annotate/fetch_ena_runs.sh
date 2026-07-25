#!/bin/bash
# Stage sequencing runs from ENA with a checksum-verified manifest.
#
# Written because ad-hoc curl loops leave no record of what was downloaded or
# whether it arrived intact. This writes a manifest first, then verifies every
# file against the ENA-published MD5, and is safe to re-run: already-verified
# files are skipped rather than re-fetched.
#
# Files are named <alias>.<run>.fastq.gz so downstream output carries the
# sample alias used in the source publication rather than a bare accession.
#
# Usage:
#   annotate/fetch_ena_runs.sh <accession> <outdir> [alias_regex]
#
#   accession    ENA/SRA project or study accession (e.g. PRJNA608699)
#   outdir       staging directory, created if absent
#   alias_regex  optional POSIX ERE matched against sample_alias to subset the
#                project (e.g. '^(MX182|SX10)$'). Omit to stage everything.
#
# Manifest columns: run_accession sample_alias read_count base_count
#                   fastq_bytes fastq_md5 fastq_ftp

set -uo pipefail

ACC="${1:?usage: fetch_ena_runs.sh <accession> <outdir> [alias_regex]}"
OUT="${2:?outdir required}"
FILTER="${3:-.}"

API="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="run_accession,sample_alias,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp"

command -v curl >/dev/null || { echo "missing dependency: curl" >&2; exit 1; }
command -v md5sum >/dev/null || { echo "missing dependency: md5sum" >&2; exit 1; }

mkdir -p "$OUT" || exit 1
cd "$OUT" || exit 1

echo "querying ENA for ${ACC} ..."
curl -sSf "${API}?accession=${ACC}&result=read_run&fields=${FIELDS}&format=tsv&limit=0" \
  | awk -F'\t' -v re="$FILTER" 'NR==1 || $2 ~ re' > manifest.tsv || {
      echo "ENA query failed for ${ACC}" >&2; exit 1; }

n=$(( $(wc -l < manifest.tsv) - 1 ))
[ "$n" -gt 0 ] || { echo "no runs matched filter '${FILTER}' in ${ACC}" >&2; exit 1; }
awk -F'\t' -v n="$n" 'NR>1{b+=$5} END{printf "%d runs, %.1f MB to stage\n", n, b/1048576}' manifest.tsv

fail=0
while IFS=$'\t' read -r run alias rc bc bytes md5 ftp; do
  [ "$run" = "run_accession" ] && continue
  [ -n "${ftp:-}" ] || { echo "SKIP  $alias ($run): no fastq_ftp published"; continue; }
  out="${alias:-$run}.${run}.fastq.gz"

  if [ -s "$out" ] && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$md5" ]; then
    echo "SKIP  $out (already verified)"
    continue
  fi

  echo "GET   $out"
  if ! curl -sSL --retry 5 --retry-delay 3 -o "$out" "https://${ftp}"; then
    echo "FAIL  download $out" >&2; fail=1; continue
  fi
  got=$(md5sum "$out" | cut -d' ' -f1)
  if [ "$got" = "$md5" ]; then
    echo "OK    $out  $got"
  else
    echo "BAD   $out expected=$md5 got=$got" >&2; fail=1
  fi
done < manifest.tsv

echo "---- staging complete, failures=${fail} ----"
exit "$fail"
