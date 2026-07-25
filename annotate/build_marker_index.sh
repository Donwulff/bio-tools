#!/usr/bin/env bash
# Build a compact marker-name -> position index from the YBrowse GFF3.
#
# Why: the YBrowse GFF3 is ~800 MB / 3.1M lines and the VCF is no better. Every
# lookup by marker *name* meant a full scan, so in practice all our tooling
# selected markers by ISOGG *label* instead. That silently drops every marker
# whose isogg_haplogroup is "unknown" / "not listed" -- which is most of what
# FTDNA's Block Tree and YFull actually use below the ISOGG-covered backbone.
#
# This index makes name-based lookup instant, so a block of SNP names copied
# from FTDNA/YFull can be tested directly.
#
# Usage: annotate/build_marker_index.sh [gff3] [out.tsv.gz]
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
gff3="${1:-$repo_root/resources/snps_hg38.gff3}"
out="${2:-$repo_root/resources/marker_index.tsv.gz}"

[[ -r "$gff3" ]] || { echo "no such GFF3: $gff3" >&2; exit 1; }

echo "indexing $gff3 -> $out" >&2

{
  printf 'name\tchrom\tpos\tanc\tder\tisogg\tycc\tyfull_node\tref\tcomment\n'
  awk -F'\t' '
    $1 ~ /^#/ { next }
    NF < 9   { next }
    {
      n=split($9, kv, ";")
      delete a
      for (i = 1; i <= n; i++) {
        eq = index(kv[i], "=")
        if (eq) a[substr(kv[i], 1, eq-1)] = substr(kv[i], eq+1)
      }
      if (a["Name"] == "") next
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        a["Name"], $1, $4,
        (a["allele_anc"] == "" ? "." : a["allele_anc"]),
        (a["allele_der"] == "" ? "." : a["allele_der"]),
        (a["isogg_haplogroup"] == "" ? "." : a["isogg_haplogroup"]),
        (a["ycc_haplogroup"] == "" ? "." : a["ycc_haplogroup"]),
        (a["yfull_node"] == "" ? "." : a["yfull_node"]),
        (a["ref"] == "" ? "." : a["ref"]),
        (a["comment"] == "" ? "." : a["comment"])
    }' "$gff3"
} | gzip -c > "$out"

echo "wrote $out ($(zcat "$out" | wc -l) rows)" >&2
