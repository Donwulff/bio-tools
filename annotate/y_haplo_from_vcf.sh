#!/bin/sh
set -eu

# Annotate chrY calls against YBrowse/ISOGG markers and report derived markers.
# Intended for haploid chrY calls, but accepts diploid 0/0 and 1/1 encodings too.

usage() {
  cat <<'EOF'
Usage:
  y_haplo_from_vcf.sh -i sample.vcf[.gz] -o out_prefix [-s SAMPLE]
                      [--ann-vcf snps_hg38.vcf.gz] [--refresh-ann]
                      [--min-gq N] [--min-dp N]
                      [--site-filter-mode any|pass|pass-or-dot|deepvariant]
                      [--site-filter-expr EXPR]
                      [--keep-tmp]

Outputs:
  <prefix>.merged.vcf.gz           # YBrowse + sample merged on chrY SNP sites
  <prefix>.called.vcf.gz           # merged VCF with HARRY/ALIEN dropped
  <prefix>.derived.vcf.gz          # sample calls derived vs INFO/AA
  <prefix>.novel_snps.vcf.gz       # sample chrY SNP ALT calls not in marker list
  <prefix>.novel_indels.vcf.gz     # sample chrY indel ALT calls not in marker list
  <prefix>.derived.tsv             # tabular derived marker list
  <prefix>.summary.txt             # deepest HG candidates + marker counts

Notes:
  - By default, marker VCF is cached at annotate/data/snps_hg38.vcf.gz.
  - --refresh-ann forces marker refresh from ybrowse.org.
  - For DeepVariant VCF, use --site-filter-mode deepvariant.
EOF
}

echo "WARN: y_haplo_from_vcf.sh is legacy. Prefer y_haplo_from_markers.py + y_path_rank.py (via run_* drivers)." >&2

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing required command: $1" >&2
    exit 1
  }
}

INPUT_VCF=""
OUT_PREFIX=""
SAMPLE=""
ANN_VCF=""
REFRESH_ANN=0
MIN_GQ=""
MIN_DP=""
KEEP_TMP=0
SITE_FILTER_MODE="any"
SITE_FILTER_EXPR=""

while [ $# -gt 0 ]; do
  case "$1" in
    -i|--input) INPUT_VCF="$2"; shift 2 ;;
    -o|--out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    -s|--sample) SAMPLE="$2"; shift 2 ;;
    --ann-vcf) ANN_VCF="$2"; shift 2 ;;
    --refresh-ann) REFRESH_ANN=1; shift ;;
    --min-gq) MIN_GQ="$2"; shift 2 ;;
    --min-dp) MIN_DP="$2"; shift 2 ;;
    --site-filter-mode) SITE_FILTER_MODE="$2"; shift 2 ;;
    --site-filter-expr) SITE_FILTER_EXPR="$2"; shift 2 ;;
    --keep-tmp) KEEP_TMP=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$INPUT_VCF" ] || { echo "ERROR: --input is required" >&2; usage; exit 1; }
[ -n "$OUT_PREFIX" ] || { echo "ERROR: --out-prefix is required" >&2; usage; exit 1; }
[ -f "$INPUT_VCF" ] || { echo "ERROR: input not found: $INPUT_VCF" >&2; exit 1; }

need_cmd bcftools
need_cmd bgzip
need_cmd tabix
need_cmd awk
need_cmd sort
need_cmd wget

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
DATA_DIR="${SCRIPT_DIR}/data"
mkdir -p "$DATA_DIR"

if [ -z "$ANN_VCF" ]; then
  ANN_VCF="${DATA_DIR}/snps_hg38.vcf.gz"
fi

ANN_TBI="${ANN_VCF}.tbi"
ANN_URL_BASE="http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz"

if [ "$REFRESH_ANN" -eq 1 ] || [ ! -f "$ANN_VCF" ] || [ ! -f "$ANN_TBI" ]; then
  echo "Refreshing marker set from YBrowse into: $ANN_VCF"
  wget -q -O "$ANN_VCF" "$ANN_URL_BASE"
  wget -q -O "$ANN_TBI" "${ANN_URL_BASE}.tbi"
fi

TMP_DIR=$(mktemp -d)
cleanup() {
  if [ "$KEEP_TMP" -eq 0 ]; then
    rm -rf "$TMP_DIR"
  else
    echo "Keeping temp files in: $TMP_DIR"
  fi
}
trap cleanup EXIT INT TERM

INPUT_GZ="$INPUT_VCF"
if gzip -t "$INPUT_VCF" >/dev/null 2>&1; then
  # Normalize gzip/bgzip input to BGZF so tabix indexing is always valid.
  INPUT_GZ="${TMP_DIR}/input.vcf.gz"
  gzip -dc "$INPUT_VCF" | bgzip -c > "$INPUT_GZ"
else
  INPUT_GZ="${TMP_DIR}/input.vcf.gz"
  bgzip -c "$INPUT_VCF" > "$INPUT_GZ"
fi

if [ ! -f "${INPUT_GZ}.tbi" ] && [ ! -f "${INPUT_GZ}.csi" ]; then
  tabix -f "$INPUT_GZ"
fi
if [ ! -f "${ANN_VCF}.tbi" ] && [ ! -f "${ANN_VCF}.csi" ]; then
  tabix -f "$ANN_VCF"
fi

if [ -z "$SAMPLE" ]; then
  SAMPLE=$(bcftools query -l "$INPUT_GZ" | head -n1)
fi
[ -n "$SAMPLE" ] || { echo "ERROR: could not determine sample name, use --sample" >&2; exit 1; }
echo "Using sample: $SAMPLE"

CHR_Y_SNPS="${TMP_DIR}/sample.chrY.snps.vcf.gz"
VIEW_EXPR=""
if [ -n "$SITE_FILTER_EXPR" ]; then
  VIEW_EXPR="$SITE_FILTER_EXPR"
else
  case "$SITE_FILTER_MODE" in
    any) VIEW_EXPR="" ;;
    pass) VIEW_EXPR='FILTER="PASS"' ;;
    pass-or-dot) VIEW_EXPR='FILTER="PASS" || FILTER="."' ;;
    deepvariant) VIEW_EXPR='FILTER="PASS" || FILTER="RefCall"' ;;
    *)
      echo "ERROR: --site-filter-mode must be one of: any, pass, pass-or-dot, deepvariant" >&2
      exit 1
      ;;
  esac
fi

if [ -n "$VIEW_EXPR" ]; then
  bcftools view -r chrY -s "$SAMPLE" -m2 -M2 -v snps -i "$VIEW_EXPR" -Oz -o "$CHR_Y_SNPS" "$INPUT_GZ"
else
  bcftools view -r chrY -s "$SAMPLE" -m2 -M2 -v snps -Oz -o "$CHR_Y_SNPS" "$INPUT_GZ"
fi
tabix -f "$CHR_Y_SNPS"

MERGED_VCF="${OUT_PREFIX}.merged.vcf.gz"
CALLED_VCF="${OUT_PREFIX}.called.vcf.gz"
DERIVED_VCF="${OUT_PREFIX}.derived.vcf.gz"
NOVEL_SNP_VCF="${OUT_PREFIX}.novel_snps.vcf.gz"
NOVEL_INDEL_VCF="${OUT_PREFIX}.novel_indels.vcf.gz"
DERIVED_TSV="${OUT_PREFIX}.derived.tsv"
SUMMARY_TXT="${OUT_PREFIX}.summary.txt"

bcftools merge -r chrY -m both -Oz -o "$MERGED_VCF" "$ANN_VCF" "$CHR_Y_SNPS"
tabix -f "$MERGED_VCF"

# Drop marker pseudo-samples and retain caller sample only.
bcftools view -s "^HARRY,ALIEN" -Oz -o "$CALLED_VCF" "$MERGED_VCF"
tabix -f "$CALLED_VCF"

DERIVED_EXPR='INFO/AA!="." && ((GT="1" && INFO/AA=REF) || (GT="1/1" && INFO/AA=REF) || (GT="0" && INFO/AA=ALT) || (GT="0/0" && INFO/AA=ALT))'
if [ -n "$MIN_GQ" ]; then
  DERIVED_EXPR="${DERIVED_EXPR} && FMT/GQ>=${MIN_GQ}"
fi
if [ -n "$MIN_DP" ]; then
  DERIVED_EXPR="${DERIVED_EXPR} && FMT/DP>=${MIN_DP}"
fi

bcftools view -i "$DERIVED_EXPR" -Oz -o "$DERIVED_VCF" "$CALLED_VCF"
tabix -f "$DERIVED_VCF"

bcftools view -i 'ID="." && TYPE="snp" && (GT="1" || GT="1/1")' -Oz -o "$NOVEL_SNP_VCF" "$CALLED_VCF"
tabix -f "$NOVEL_SNP_VCF"

bcftools view -i 'ID="." && TYPE="indel" && (GT="1" || GT="1/1")' -Oz -o "$NOVEL_INDEL_VCF" "$CALLED_VCF"
tabix -f "$NOVEL_INDEL_VCF"

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AA\t%INFO/HG\t%INFO/ISOGG\t[%GT\t%GQ\t%DP]\n' "$DERIVED_VCF" \
  > "$DERIVED_TSV"

{
  CALLED_NONMISS=$(bcftools view -i 'GT!="mis"' -H "$CALLED_VCF" | wc -l)
  CALLED_TOTAL=$(bcftools view -H "$CALLED_VCF" | wc -l)

  echo "Input: $INPUT_VCF"
  echo "Sample: $SAMPLE"
  echo "Markers: $ANN_VCF"
  echo "MinGQ: ${MIN_GQ:-none}"
  echo "MinDP: ${MIN_DP:-none}"
  echo "SiteFilterMode: $SITE_FILTER_MODE"
  echo "SiteFilterExpr: ${SITE_FILTER_EXPR:-none}"
  echo
  echo "Counts:"
  printf "  merged:         %s\n" "$(bcftools view -H "$MERGED_VCF" | wc -l)"
  printf "  called_total:   %s\n" "$CALLED_TOTAL"
  printf "  called_nonmiss: %s\n" "$CALLED_NONMISS"
  printf "  derived:        %s\n" "$(bcftools view -H "$DERIVED_VCF" | wc -l)"
  printf "  novel_snps:     %s\n" "$(bcftools view -H "$NOVEL_SNP_VCF" | wc -l)"
  printf "  novel_indels:   %s\n" "$(bcftools view -H "$NOVEL_INDEL_VCF" | wc -l)"
  if [ "$CALLED_NONMISS" -eq 0 ]; then
    echo
    echo "WARNING:"
    echo "  No non-missing genotypes at marker sites."
    echo "  Input is likely variant-only VCF. Use a VCF/gVCF with ref-calls at marker positions,"
    echo "  or genotype marker loci directly from BAM."
  fi
  echo
  echo "Deepest HG candidates (max chain depth by '>'):"
  awk -F'\t' '
    {
      hg=$7
      gsub(/"/, "", hg)
      if (hg == "" || hg == ".") next
      tmp=hg
      depth=gsub(/>/, ">", tmp)+1
      if (depth > max_depth) {
        max_depth=depth
        delete best
        best[hg]=1
      } else if (depth == max_depth) {
        best[hg]=1
      }
    }
    END {
      if (max_depth == 0) {
        print "  (none)"
      } else {
        printf "  depth=%d\n", max_depth
        for (h in best) print "  " h
      }
    }
  ' "$DERIVED_TSV"
  echo
  echo "Top ISOGG labels among derived markers:"
  awk -F'\t' '
    {
      g=$8
      gsub(/"/, "", g)
      if (g == "" || g == ".") g="(missing)"
      c[g]++
    }
    END {
      for (k in c) printf "%d\t%s\n", c[k], k
    }
  ' "$DERIVED_TSV" | sort -nr | head -n 20 | awk '{print "  " $0}'
} > "$SUMMARY_TXT"

echo "Done. Summary: $SUMMARY_TXT"
