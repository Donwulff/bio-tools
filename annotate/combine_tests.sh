#!/bin/bash
set -v -x

# In case it's necessary to set paths etc. these aren't modified but vcf conversions are added
TSV=genome_v3_Full
GFG=GFG_annotated_unphased_genotypes
ANC=AncestryDNA

HS37D5=/mnt/GenomicData/hs37d5/hs37d5.fa
WORK=work

mkdir -p ${WORK}

convert_23andMe()
{
  # WARNING: This is very incomplete, some 23andMe SNV's are on wrong orientation!
  # bcftools convert --tsv2vcf only drops insertions and deletions:
  # join --nocheck-order -j1 -v1  <(grep -v "^#" genome_v3_Full.txt | gawk '{ print $2":"$3,$1,$4 }') <(zgrep -v "^#" genome_v3_Full.norm.vcf.gz \
  #   | gawk '{ print $1":"$2 }') | egrep -v "(I|D)"
  bcftools convert -f ${HS37D5} --tsv2vcf ${TSV}.txt -c ID,CHROM,POS,AA --samples PARTICIPANT 2> ${TSV}.log \
    | bcftools norm -f ${HS37D5} -cws 2> ${TSV}.norm.log \
    | bgzip -c > ${TSV}.norm.vcf.gz
  tabix -f ${TSV}.norm.vcf.gz
}

convert_Genes4Good()
{
  # bcftols norm will write zero length contigs if given file, so pipe it, and also filter out contigs just in case
  # It also marks deletions with - as alt, which makes most downstream analysis choke.
  # Running norm after blanking out the -'s causes them to turn into SNV's, so avoid that.
  zcat ${GFG}.vcf.gz \
    | bcftools norm -f ${HS37D5} -cws 2> ${GFG}.norm.log \
    | gawk '{ OFS="\t"; if ($5=="-") $5="."; print }' \
    | grep -v "^##contig" \
    | bgzip -c > ${GFG}.norm.vcf.gz
  tabix -f ${GFG}.norm.vcf.gz
}

convert_AncestryDNA()
{
  # Fix the chromosome ploidy for male sample
  cat ${ANC}.txt | gawk '!/^#/{if(0!=$4){if(23==$2||24==$2||25==$2){if(substr($4,0,1)!=substr($5,0,1))$4="-";$5=""};print $1,$2,$3,$4$5}}' > ${ANC}.pld
  # Fix the chromosome names
  cat ${ANC}.pld | gawk '!/^#/{if(0!=$4){if(23==$2||25==$2)$2="X";if(24==$2)$2="Y";if(26==$2)$2="MT";print $1,$2,$3,$4$5}}' > ${ANC}.chr
  # unfortunately AncestryDNA has PAR1 after Y, so we need to do some magic.
  cat ${ANC}.chr | gawk '{if($2!="X"&&$2!="Y"&&$2!="MT")print}' > ${ANC}.new
  cat ${ANC}.chr | gawk '{if($2=="X")print}' | sort -k1,1V -k2,2n >> ${ANC}.new
  cat ${ANC}.chr | gawk '{if($2=="Y"||$2=="MT")print}' >> ${ANC}.new
  # bcftools norm: NON_ACGTN appears to be deleted with - alt or ref. -cws's s doesn't fix them. Leave for now, if they're indels
  bcftools convert -f ${HS37D5} --tsv2vcf ${ANC}.new -c ID,CHROM,POS,AA --samples PARTICIPANT 2> ${ANC}.log \
    | bcftools norm -f ${HS37D5} -cws 2> ${ANC}.norm.log \
    | bgzip -c > ${ANC}.norm.vcf.gz
  tabix -f ${ANC}.norm.vcf.gz
}

convert_23andMe
convert_Genes4Good
convert_AncestryDNA

# Put 23andMe conversion first, since it has correct contigs; -D leaves in all duplicates, no option to delete
bcftools concat -a -D ${TSV}.norm.vcf.gz ${GFG}.norm.vcf.gz ${ANC}.norm.vcf.gz | bgzip -c > ${WORK}/merged.vcf.gz

# NB. if tested variants are different and thus the results differ, they're removed.
# These should probably be accepted for the the different alt, but it's tricky (What if both variants are 1/0?)
# Yugh, diff-site-discordance seems to take first matching location only, so we need to do this hard way
> ${WORK}/discordant.sites
vcftools --gzvcf ${TSV}.norm.vcf.gz --gzdiff ${GFG}.norm.vcf.gz --diff-site-discordance --out ${WORK}/merged
grep -v -- "-nan" ${WORK}/merged.diff.sites  | grep -v "0$" >> ${WORK}/discordant.sites
vcftools --gzvcf ${TSV}.norm.vcf.gz --gzdiff ${ANC}.norm.vcf.gz --diff-site-discordance --out ${WORK}/merged
grep -v -- "-nan" ${WORK}/merged.diff.sites  | grep -v "0$" >> ${WORK}/discordant.sites
vcftools --gzvcf ${GFG}.norm.vcf.gz --gzdiff ${ANC}.norm.vcf.gz --diff-site-discordance --out ${WORK}/merged
grep -v -- "-nan" ${WORK}/merged.diff.sites  | grep -v "0$" >> ${WORK}/discordant.sites
# Use --remove-indels if needed, but imputation should just ignore it, variant browsers and annotators may be able to use them.
vcftools --gzvcf ${WORK}/merged.vcf.gz --exclude-positions ${WORK}/discordant.sites --out ${WORK}/filtered --recode --stdout | bgzip -c > filtered.vcf.gz
tabix -f filtered.vcf.gz

# Some tools require split by chromosome/contig
#mkdir -p chrom
#tabix -l filtered.vcf.gz | xargs -IX bash -c "tabix -h filtered.vcf.gz X | bgzip -c > chrom/X.vcf.gz"
