#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  sanitize_marker_vcf.sh --input in.vcf.gz --output out.vcf.gz

Sanitizes marker VCF alleles for strict tools (e.g. LiftoverVcf):
- removes ALT alleles equal to REF
- removes duplicate ALT alleles in the same record
- drops records with no remaining ALT alleles
EOF
}

INPUT=""
OUTPUT=""

while [ $# -gt 0 ]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$INPUT" ] || { echo "ERROR: --input is required" >&2; exit 1; }
[ -n "$OUTPUT" ] || { echo "ERROR: --output is required" >&2; exit 1; }
[ -f "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools missing" >&2; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip missing" >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix missing" >&2; exit 1; }

TMP="${OUTPUT}.tmp"
mkdir -p "$(dirname "$OUTPUT")"

bcftools view -Ov "$INPUT" \
  | awk 'BEGIN{OFS="\t"; fixed=0; dropped=0}
         /^#/ { print; next }
         {
           ref=$4
           n=split($5, a, ",")
           m=0
           delete seen
           delete b
           changed=0
           for (i=1; i<=n; i++) {
             alt=a[i]
             if (alt=="." || alt==ref) { changed=1; continue }
             if (seen[alt]) { changed=1; continue }
             seen[alt]=1
             b[++m]=alt
           }
           if (m==0) { dropped++; next }
           if (changed || m!=n) { fixed++ }
           $5=b[1]
           for (i=2; i<=m; i++) $5=$5","b[i]
           print
         }
         END{
           printf("Sanitized marker VCF: fixed=%d dropped=%d\n", fixed, dropped) > "/dev/stderr"
         }' \
  | bgzip -c > "$TMP"

mv -f "$TMP" "$OUTPUT"
tabix -f "$OUTPUT"
