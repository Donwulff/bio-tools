#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  samtools_if.sh <subcmd> [samtools-args...] <in.bam>
  samtools_if.sh <subcmd> [samtools-args...] -- <in.bam>

Behavior:
  - Writes output to <in.bam>.<subcmd> by default.
  - Skips work if output exists and is newer than <in.bam>.
  - Use -o/--out to override output path.
  - Use -f/--force to re-run even if up to date.

Examples:
  samtools_if.sh flagstat sample.bam
  samtools_if.sh stats -r chr1 -- sample.bam
  samtools_if.sh idxstats -o sample.idxstats sample.bam
USAGE
}

if [[ $# -lt 2 ]]; then
  usage
  exit 2
fi

subcmd="$1"
shift

out=""
force=0
args=()
infile=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--out)
      out="$2"
      shift 2
      ;;
    -f|--force)
      force=1
      shift
      ;;
    --)
      shift
      infile="${1:-}"
      shift || true
      args+=("$@")
      break
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

if [[ -z "$infile" ]]; then
  if [[ ${#args[@]} -eq 0 ]]; then
    usage
    exit 2
  fi
  infile="${args[-1]}"
  args=("${args[@]:0:${#args[@]}-1}")
fi

if [[ -z "$infile" || ! -f "$infile" ]]; then
  echo "ERROR: input file not found: $infile" >&2
  exit 2
fi

if [[ -z "$out" ]]; then
  out="${infile}.${subcmd}"
fi

if [[ $force -eq 0 && -s "$out" && "$out" -nt "$infile" ]]; then
  echo "Up to date: $out" >&2
  exit 0
fi

samtools "$subcmd" "${args[@]}" "$infile" | tee "$out"
