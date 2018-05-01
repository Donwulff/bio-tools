#!/bin/bash
set -v

TSV=genome_v3_Full
GFG=GFG_annotated_unphased_genotypes

HS37D5=/mnt/GenomicData/hs37d5/hs37d5.fa

convert_23andMe()
{
  # WARNING: This is very incomplete, some 23andMe SNV's are on wrong orientation!
  bcftools convert -f ${HS37D5} --tsv2vcf ${TSV}.txt -c ID,CHROM,POS,AA --samples PARTICIPANT 2> ${TSV}.log \
    | bcftools norm -f ${HS37D5} -cws 2> ${TSV}.norm.log \
    | bgzip -c > ${TSV}.norm.vcf.gz
  tabix -f ${TSV}.norm.vcf.gz
}

convert_genes4good()
{
  # bcftols norm will write zero length contigs if given file, so pipe it, and also filter out contigs just in case
  # It also marks deletions with - as alt, which makes most downstream analysis choke.
  # Running norm after filtering out the -'s causes them to turn into SNP's, so avoid that.
  zcat ${GFG}.vcf.gz \
    | bcftools norm -f ${HS37D5} -cws 2> ${GFG}.norm.log \
    | gawk '{ OFS="\t"; if ($5=="-") $5="."; print }' \
    | grep -v "^##contig" \
    | bgzip -c > ${GFG}.norm.vcf.gz
  tabix -f ${GFG}.norm.vcf.gz
}

convert_23andMe
convert_genes4good

# Put 23andMe conversion first, since it has correct contigs; delete duplicates.
bcftools concat -a -D ${TSV}.norm.vcf.gz ${GFG}.norm.vcf.gz | bgzip -c > merged.vcf.gz

# NB. if tested variants are different and thus the results differ, they're removed.
# These should probably be accepted for the the different alt, but it's tricky (What if both variants are 1/0?)
vcftools --gzvcf ${TSV}.norm.vcf.gz --gzdiff ${GFG}.norm.vcf.gz --diff-site-discordance --out merged
grep -v -- "-nan" merged.diff.sites  | grep -v "0$" > discordant.sites
vcftools --gzvcf merged.vcf.gz --exclude-positions discordant.sites --out filtered --recode --stdout | bgzip -c > filtered.vcf.gz
