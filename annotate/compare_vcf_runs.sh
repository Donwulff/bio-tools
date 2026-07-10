#!/bin/sh
set -eu

usage() {
  cat <<'EOF'
Usage:
  compare_vcf_runs.sh -a runA.vcf[.gz] -b runB.vcf[.gz] [options]

Compares two single-sample VCF callsets from the same sample/reference.
Designed for DeepVariant branch comparisons (e.g. small-model on/off).

Options:
  -a, --vcf-a FILE          First VCF (required)
  -b, --vcf-b FILE          Second VCF (required)
  --label-a NAME            Label for first VCF (default: a)
  --label-b NAME            Label for second VCF (default: b)
  -o, --out-dir DIR         Output directory (default: ./vcf_compare_YYYYmmdd_HHMMSS)
  --regions-file FILE       Optional TSV: CHROM<TAB>POS<TAB>LABEL(optional)
  --keep-work               Keep normalized working files
  -h, --help                Show this help

Outputs:
  summary.tsv                     # all + PASS-level site/GT concordance summary
  all.only_<label>.vcf.gz         # private sites for each run
  all.shared_<label>.vcf.gz       # shared sites represented from each run
  all.gt_compare.tsv              # shared-site GT concordance table
  pass.*                          # same metrics restricted to FILTER=PASS
  regions_compare.tsv             # optional per-site comparison for --regions-file
EOF
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

VCF_A=""
VCF_B=""
LABEL_A="a"
LABEL_B="b"
OUT_DIR="./vcf_compare_$(date +%Y%m%d_%H%M%S)"
REGIONS_FILE=""
KEEP_WORK=0

while [ $# -gt 0 ]; do
  case "$1" in
    -a|--vcf-a) VCF_A="$2"; shift 2 ;;
    -b|--vcf-b) VCF_B="$2"; shift 2 ;;
    --label-a) LABEL_A="$2"; shift 2 ;;
    --label-b) LABEL_B="$2"; shift 2 ;;
    -o|--out-dir) OUT_DIR="$2"; shift 2 ;;
    --regions-file) REGIONS_FILE="$2"; shift 2 ;;
    --keep-work) KEEP_WORK=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$VCF_A" ] || { echo "ERROR: --vcf-a required" >&2; usage; exit 1; }
[ -n "$VCF_B" ] || { echo "ERROR: --vcf-b required" >&2; usage; exit 1; }
[ -f "$VCF_A" ] || { echo "ERROR: missing file: $VCF_A" >&2; exit 1; }
[ -f "$VCF_B" ] || { echo "ERROR: missing file: $VCF_B" >&2; exit 1; }
[ -z "$REGIONS_FILE" ] || [ -f "$REGIONS_FILE" ] || { echo "ERROR: missing file: $REGIONS_FILE" >&2; exit 1; }

need_cmd bcftools
need_cmd bgzip
need_cmd tabix
need_cmd awk
need_cmd wc
need_cmd paste

mkdir -p "$OUT_DIR"
WORK_DIR="$OUT_DIR/work"
mkdir -p "$WORK_DIR"

cleanup() {
  if [ "$KEEP_WORK" -eq 0 ]; then
    rm -rf "$WORK_DIR"
  fi
}
trap cleanup EXIT INT TERM

normalize_vcf() {
  in_vcf="$1"
  out_vcfgz="$2"
  # Handles plain VCF, bgzip VCF, and gzip-compressed VCF even without .gz suffix.
  bcftools view -Oz -o "$out_vcfgz" "$in_vcf"
  bcftools index -f -t "$out_vcfgz"
}

compare_set() {
  set_name="$1"
  vcf_a="$2"
  vcf_b="$3"

  only_a="$OUT_DIR/${set_name}.only_${LABEL_A}.vcf.gz"
  only_b="$OUT_DIR/${set_name}.only_${LABEL_B}.vcf.gz"
  shared_a="$OUT_DIR/${set_name}.shared_${LABEL_A}.vcf.gz"
  shared_b="$OUT_DIR/${set_name}.shared_${LABEL_B}.vcf.gz"
  gt_a="$WORK_DIR/${set_name}.gt_${LABEL_A}.tsv"
  gt_b="$WORK_DIR/${set_name}.gt_${LABEL_B}.tsv"
  gt_cmp="$OUT_DIR/${set_name}.gt_compare.tsv"

  bcftools isec -n=1 -w1 -Oz -o "$only_a" "$vcf_a" "$vcf_b"
  bcftools isec -n=1 -w2 -Oz -o "$only_b" "$vcf_a" "$vcf_b"
  bcftools isec -n=2 -w1 -Oz -o "$shared_a" "$vcf_a" "$vcf_b"
  bcftools isec -n=2 -w2 -Oz -o "$shared_b" "$vcf_a" "$vcf_b"

  tabix -f "$only_a" >/dev/null
  tabix -f "$only_b" >/dev/null
  tabix -f "$shared_a" >/dev/null
  tabix -f "$shared_b" >/dev/null

  total_a=$(bcftools view -H "$vcf_a" | wc -l | tr -d ' ')
  total_b=$(bcftools view -H "$vcf_b" | wc -l | tr -d ' ')
  n_only_a=$(bcftools view -H "$only_a" | wc -l | tr -d ' ')
  n_only_b=$(bcftools view -H "$only_b" | wc -l | tr -d ' ')
  n_shared=$(bcftools view -H "$shared_a" | wc -l | tr -d ' ')

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$shared_a" > "$gt_a"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$shared_b" > "$gt_b"

  paste "$gt_a" "$gt_b" > "$gt_cmp"

  gt_same=$(awk 'BEGIN{FS="\t"} $1==$6 && $2==$7 && $3==$8 && $4==$9 && $5==$10 {c++} END{print c+0}' "$gt_cmp")
  gt_diff=$(awk 'BEGIN{FS="\t"} $1==$6 && $2==$7 && $3==$8 && $4==$9 && $5!=$10 {c++} END{print c+0}' "$gt_cmp")
  gt_diff_pct=$(awk -v d="$gt_diff" -v s="$n_shared" 'BEGIN{if(s>0) printf "%.6f", 100.0*d/s; else printf "0.000000"}')

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$set_name" "$LABEL_A" "$LABEL_B" "$total_a" "$total_b" "$n_shared" "$n_only_a" "$n_only_b" "$gt_same" "$gt_diff" \
    >> "$WORK_DIR/summary_rows.tsv"
  printf '%s\t%s\n' "${set_name}_shared_gt_diff_pct" "$gt_diff_pct" >> "$WORK_DIR/summary_kv.tsv"
}

query_site_row() {
  vcf="$1"
  chrom="$2"
  pos="$3"
  row=$(bcftools query -r "${chrom}:${pos}-${pos}" -f '[%GT\t%GQ\t%DP]\t%QUAL\t%FILTER\t%REF\t%ALT\n' "$vcf" | head -n1 || true)
  if [ -z "$row" ]; then
    printf '.\t.\t.\t.\t.\t.\t.'
  else
    printf '%s' "$row"
  fi
}

NORM_A="$WORK_DIR/${LABEL_A}.norm.vcf.gz"
NORM_B="$WORK_DIR/${LABEL_B}.norm.vcf.gz"
normalize_vcf "$VCF_A" "$NORM_A"
normalize_vcf "$VCF_B" "$NORM_B"

printf 'set\tlabel_a\tlabel_b\ttotal_a\ttotal_b\tshared\tonly_a\tonly_b\tshared_gt_same\tshared_gt_diff\n' > "$OUT_DIR/summary.tsv"
: > "$WORK_DIR/summary_rows.tsv"
: > "$WORK_DIR/summary_kv.tsv"

compare_set "all" "$NORM_A" "$NORM_B"

PASS_A="$WORK_DIR/${LABEL_A}.pass.vcf.gz"
PASS_B="$WORK_DIR/${LABEL_B}.pass.vcf.gz"
bcftools view -f PASS -Oz -o "$PASS_A" "$NORM_A"
bcftools view -f PASS -Oz -o "$PASS_B" "$NORM_B"
tabix -f "$PASS_A" >/dev/null
tabix -f "$PASS_B" >/dev/null

compare_set "pass" "$PASS_A" "$PASS_B"

cat "$WORK_DIR/summary_rows.tsv" >> "$OUT_DIR/summary.tsv"

{
  printf 'metric\tvalue\n'
  cat "$WORK_DIR/summary_kv.tsv"
} > "$OUT_DIR/summary_metrics.tsv"

if [ -n "$REGIONS_FILE" ]; then
  out_regions="$OUT_DIR/regions_compare.tsv"
  printf 'label\tchrom\tpos\t%s_gt\t%s_gq\t%s_dp\t%s_qual\t%s_filter\t%s_ref\t%s_alt\t%s_gt\t%s_gq\t%s_dp\t%s_qual\t%s_filter\t%s_ref\t%s_alt\n' \
    "$LABEL_A" "$LABEL_A" "$LABEL_A" "$LABEL_A" "$LABEL_A" "$LABEL_A" "$LABEL_A" \
    "$LABEL_B" "$LABEL_B" "$LABEL_B" "$LABEL_B" "$LABEL_B" "$LABEL_B" "$LABEL_B" \
    > "$out_regions"

  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=2 {label=(NF>=3 && $3!="" ? $3 : $1":"$2); print $1,$2,label}' "$REGIONS_FILE" \
    | while IFS="$(printf '\t')" read -r chrom pos label; do
        ra=$(query_site_row "$NORM_A" "$chrom" "$pos")
        rb=$(query_site_row "$NORM_B" "$chrom" "$pos")
        printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$chrom" "$pos" "$ra" "$rb" >> "$out_regions"
      done
fi

echo "Done. Output dir: $OUT_DIR"
