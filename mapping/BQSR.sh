#!/bin/sh

SAMPLE=sample1-fastp-hg38.srt.bam
DATA=/mnt/GenomicData
CORES=`nproc`

# Genoe Analysis ToolKit at Broad Institute maintains a resource bundle including validated indels in human genome https://software.broadinstitute.org/gatk/download/bundle
wget -nc ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
if [ -e ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.nohla.vcf.gz.tbi ];
then
  # Unfortunately the header contigs include empty HLA contigs which depend on HLA release version, so we need to cut them out or GATK throws a fit
  tabix -H ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | grep -v "^##contig=<ID=HLA-" > Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.head
  bcftools reheader -h Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.head ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -o ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.nohla.vcf.gz
  tabix ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.nohla.vcf.gz
fi

# dbSNP database snapshot, version 151, GATK contig names already provided by dbSNP at National Center for Biotechnology Information NCBI, National Institute of Health https://www.ncbi.nlm.nih.gov/projects/SNP/
wget -nc ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz -O ${DATA}/All_20180418_GRCh38p7_GATK.vcf.gz
tabix ${DATA}/All_20180418_GRCh38p7_GATK.vcf.gz

if [ -e ${SAMPLE}.recal ];
then
  tabix -l ${DATA}/All_20180418_GRCh38p7_GATK.vcf.gz \
    | xargs -IZ -P${CORES} gatk-4.0.4.0/gatk --java-options -Xms4G BaseRecalibrator -R /mnt/GenomicData/GRCh38/hg38.fa \
      --known-sites ${DATA}\Mills_and_1000G_gold_standard.indels.hg38.nohla.vcf.gz \
      --known-sites ${DATA}\All_20180418_GRCh38p7_GATK.vcf.gz \
      -I ${SAMPLE} -O ${SAMPLE}.Z.recal -L Z

  # Tests show that the male X and Y chromosome covariates differ from the rest of the genome, perhaps due to being phased. (Test this hypothesis on female X later)
  # Because we have no reason to expect the sequencing error profile to differ by chromosome, process autosomal chromosomes only.
  ls ${SAMPLE}.chr[0-9]*.recal > chr.list
  gatk-4.0.4.0/gatk GatherBQSRReports -I chr.list -O ${SAMPLE}.recal
  gatk-4.0.4.0/gatk AnalyzeCovariates --bqsr ${SAMPLE}.recal --before ${SAMPLE}.chrX.recal --after ${SAMPLE}.chrY.recal --plots ${SAMPLE}.XY.pdf
fi
