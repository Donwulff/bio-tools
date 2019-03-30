#!/bin/sh

SAMPLE=${1:-sample1.srt.bam}
DATA=/mnt/GenomicData
CORES=`nproc`
GATK=4.0.5.2

if [ ! -e gatk-${GATK} ];
then
  wget -nc https://github.com/broadinstitute/gatk/releases/download/${GATK}/gatk-${GATK}.zip
  unzip gatk-${GATK}.zip
fi

set -x

# GRCh38 based with chr* chromosomes
REF=${DATA}/GRCh38/hg38.fa
DBSNP=${DATA}/All_20180418_GRCh38p7_GATK.vcf.gz
INDEL1=${DATA}/Mills_and_1000G_gold_standard.indels.hg38.noHLA.vcf.gz
INDEL2=${DATA}/Homo_sapiens_assembly38.known_indels.noHLA.vcf.gz

# hg37 based with numerical chromosomes
#REF=${DATA}/bwakit/hs37d5.fa
#DBSNP=${DATA}/All_20180423_GRCh37p13.vcf.gz
#INDEL1=${DATA}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
#INDEL2=${DATA}/1000G_phase1.indels.b37.vcf.gz

if [ ! -e ${SAMPLE%%.bam}.bai ];
then
  samtools index -@${CORES} ${SAMPLE}
fi

# Genome Analysis ToolKit at Broad Institute maintains a resource bundle including validated indels in human genome https://software.broadinstitute.org/gatk/download/bundle
wget -nc ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
if [ ! -e ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.noHLA.vcf.gz.tbi ];
then
  # Unfortunately the header contigs include empty HLA contigs which depend on HLA release version, so we need to cut them out or GATK throws a fit
  tabix -H ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | grep -v "^##contig=<ID=HLA-" > ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.noHLA.vcf.head
  bcftools reheader -h ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.noHLA.vcf.head ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                    -o ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.noHLA.vcf.gz
  tabix ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.noHLA.vcf.gz
fi
wget -nc ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -O ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
tabix ${DATA}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# dbSNP database snapshot, version 151, GATK contig names already provided by dbSNP at National Center for Biotechnology Information NCBI, National Institute of Health https://www.ncbi.nlm.nih.gov/projects/SNP/
wget -nc ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz -O ${DATA}/All_20180418_GRCh38p7_GATK.vcf.gz
tabix ${DATA}/All_20180418_GRCh38p7_GATK.vcf.gz
wget -nc ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz -O ${DATA}/All_20180423_GRCh37p13.vcf.gz
tabix ${DATA}/All_20180423_GRCh37p13.vcf.gz

# Equivalent doesn't exist for GRCh38, however the new GATK 4.0 Best Practices scripts are using known_indels file.
#wget -nc ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz -O ${DATA}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
#tabix ${DATA}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -nc ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz -O ${DATA}/Homo_sapiens_assembly38.known_indels.vcf.gz
if [ ! -e ${DATA}/Homo_sapiens_assembly38.known_indels.noHLA.vcf.gz ];
then
  # Unfortunately the header contigs include empty HLA contigs which depend on HLA release version, so we need to cut them out or GATK throws a fit
  tabix -H ${DATA}/Homo_sapiens_assembly38.known_indels.vcf.gz | grep -v "^##contig=<ID=HLA-" > ${DATA}/Homo_sapiens_assembly38.known_indels.noHLA.vcf.head
  bcftools reheader -h ${DATA}/Homo_sapiens_assembly38.known_indels.noHLA.vcf.head ${DATA}/Homo_sapiens_assembly38.known_indels.vcf.gz \
                    -o ${DATA}/Homo_sapiens_assembly38.known_indels.noHLA.vcf.gz
  tabix ${DATA}/Homo_sapiens_assembly38.known_indels.noHLA.vcf.gz
fi
wget -nc ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz -O ${DATA}/1000G_phase1.indels.b37.vcf.gz
tabix ${DATA}/1000G_phase1.indels.b37.vcf.gz

if [ ! -e ${SAMPLE}.recal ];
then
  tabix -l ${DBSNP} \
    | xargs -IZ -P${CORES} gatk-${GATK}/gatk --java-options -Xms4G BaseRecalibrator -R ${REF} \
      --known-sites ${DBSNP} \
      --known-sites ${INDEL1} \
      --known-sites ${INDEL2} \
      -I ${SAMPLE} -O ${SAMPLE}.Z.recal -L Z

  # Tests show that the male X and Y chromosome covariates differ from the rest of the genome, perhaps due to being phased. (Test this hypothesis on female X later)
  # Because we have no reason to expect the sequencing error profile to differ by chromosome, process autosomal chromosomes only.
  ls ${SAMPLE}.*[0-9].recal > chr.list
  gatk-${GATK}/gatk GatherBQSRReports -I chr.list -O ${SAMPLE}.recal
  gatk-${GATK}/gatk AnalyzeCovariates --bqsr ${SAMPLE}.recal --before ${SAMPLE}.*X.recal --after ${SAMPLE}.*Y.recal --plots ${SAMPLE}.XY.pdf
fi

if [ ! -e ${SAMPLE%%.srt.bam}.bqsr.bam ];
then
  gatk-${GATK}/gatk ApplyBQSR -R ${REF} --bqsr ${SAMPLE}.recal -I ${SAMPLE} -O ${SAMPLE%%.srt.bam}.bqsr.bam
fi

CHR=`tabix -l ${DBSNP} | grep 20`
gatk-${GATK}/gatk --java-options -Xms4G BaseRecalibrator -R ${REF} \
  --known-sites ${DBSNP} --known-sites ${INDEL1} --known-sites ${INDEL2} -I ${SAMPLE%%.srt.bam}.bqsr.bam -O ${SAMPLE}.${CHR}.after -L ${CHR}
gatk-${GATK}/gatk AnalyzeCovariates --bqsr ${SAMPLE}.recal --before ${SAMPLE}.${CHR}.recal --after ${SAMPLE}.${CHR}.after --plots ${SAMPLE}.pdf
