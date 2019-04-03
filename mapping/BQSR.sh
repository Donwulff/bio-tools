#!/bin/sh

SAMPLE=${1:-sample1.srt.bam}
DATA=/mnt/GenomicData
CORES=`nproc`
GATK=4.1.1.0

# Requires wget gawk tabix samtools bcftools python java

if [ ! -e gatk-${GATK} ];
then
  wget -nc https://github.com/broadinstitute/gatk/releases/download/${GATK}/gatk-${GATK}.zip
  unzip gatk-${GATK}.zip
fi

set -x

# GRCh38 based with chr* chromosomes
REF=${DATA}/GRCh38/hg38.fa
DBSNP=${DATA}/GCF_000001405.38.GATK.bgz
INDEL1=${DATA}/Mills_and_1000G_gold_standard.indels.hg38.noHLA.vcf.gz
INDEL2=${DATA}/Homo_sapiens_assembly38.known_indels.noHLA.vcf.gz
YBROWSE=${DATA}/snps_hg38_GATK.vcf.gz

# hg37 based with numerical chromosomes - GRCh37 is becoming obsolete, this may need some tweaking.
#REF=${DATA}/hs37d5/hs37d5.fa
#DBSNP=${DATA}/All_20180423_GRCh37p13.vcf.gz
#INDEL1=${DATA}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
#INDEL2=${DATA}/1000G_phase1.indels.b37.vcf.gz
#YBROWSE=${DATA}/snps_hg19_GATK.vcf.gz

# If bio-tools.cfg exists where the script is run from, over-ride settings from it.
if [ -e ./bio-tools.cfg ];
then
  . ./bio-tools.cfg
fi

if [ ! -e ${SAMPLE}.bai ];
then
  samtools index -@${CORES} ${SAMPLE}
fi

# GRCh37 dbSNP database snapshot, version 151, GATK contig names already provided by dbSNP at National Center for Biotechnology Information NCBI, National Institute of Health https://www.ncbi.nlm.nih.gov/projects/SNP/
wget -nc ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz -O ${DATA}/All_20180423_GRCh37p13.vcf.gz
tabix ${DATA}/All_20180423_GRCh37p13.vcf.gz

# GRCh38 dbSNP database snapshot 152, National Center for Biotechnology Information NCBI, National Institute of Health https://www.ncbi.nlm.nih.gov/projects/SNP/
# dbSNP 152 dump is 14 gigabytes, I need to devise a convention for handling files in this script.
# Right now the download goes to working directory and GATK prepared version into defined path!
if [ ! -e ${DBSNP}.tbi ];
then
  # GRCh37, now getting obsolete.
  #wget ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.bgz

  wget -nc ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.bgz
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_assembly_report.txt

  # File for next build; these will be hard to find by hand.
  #wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt

  gawk 'BEGIN { FS="\t" } !/^#/ { print $7,$10  }' GCF_000001405.38_GRCh38.p12_assembly_report.txt > NCBI-to-UCSC-GRCh38.p12.map
  time bcftools annotate --rename-chrs NCBI-to-UCSC-GRCh38.p12.map GCF_000001405.38.bgz | gawk '/^#/ && !/^##contig=/ { print } !/^#/ { if( $1!="na" ) print }' \
    | bgzip -@${CORES} -l9 -c > ${DATA}/GCF_000001405.38.GATK.bgz
  tabix ${DATA}/GCF_000001405.38.GATK.bgz
fi

# YBrowse Y-chromosome SNP list; this updates frequently, delete and redownload if needed.
# We need to filter out incorrect contig-line and variants that match reference.
# Only GRCh38 is updated, we'll need to liftover if we ever want to use new SNP's on GRCh37
wget -nc http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz
wget -nc http://ybrowse.org/gbrowse2/gff/snps_hg19.vcf.gz
if [ ! -e ${YBROWSE}.tbi ];
then
  zgrep -v "contig" snps_hg38.vcf.gz | bcftools norm -cs -f ${REF} | bgzip -c > ${DATA}/snps_hg38_GATK.vcf.gz
  tabix -f ${DATA}/snps_hg38_GATK.vcf.gz
  zgrep -v "contig" snps_hg19.vcf.gz | bcftools norm -cs -f ${REF} | bgzip -c > ${DATA}/snps_hg19_GATK.vcf.gz
  tabix -f ${DATA}/snps_hg19_GATK.vcf.gz
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

# Equivalent to this doesn't exist for GRCh38, however the new GATK 4.0 Best Practices scripts are using known_indels file.
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
      --known-sites ${YBROWSE} \
      -I ${SAMPLE} -O ${SAMPLE}.Z.recal -L Z

  # Tests show that the male X and Y chromosome covariates differ from the rest of the genome, perhaps due to being phased. (Test this hypothesis on female X later)
  # Now with YBrowse SNP's added, just use whole primary assembly covariates by default.
  ls ${SAMPLE}.*[0-9].recal > chr.list
  ls ${SAMPLE}.chrX.recal ${SAMPLE}.chrY.recal >> chr.list
  gatk-${GATK}/gatk GatherBQSRReports -I chr.list -O ${SAMPLE}.recal

  # Generate a report showing difference between X and Y chromosome BQSR covariates for checking.
  gatk-${GATK}/gatk AnalyzeCovariates --bqsr ${SAMPLE}.recal --before ${SAMPLE}.*X.recal --after ${SAMPLE}.*Y.recal --plots ${SAMPLE}.YvsX.pdf
fi

if [ ! -e ${SAMPLE%%.srt.bam}.bqsr.bam ];
then
  gatk-${GATK}/gatk ApplyBQSR -R ${REF} --bqsr ${SAMPLE}.recal -I ${SAMPLE} -O ${SAMPLE%%.srt.bam}.bqsr.bam
fi

CHR=`tabix -l ${DBSNP} | grep -m1 '20$'`
gatk-${GATK}/gatk --java-options -Xms4G BaseRecalibrator -R ${REF} \
  --known-sites ${DBSNP} --known-sites ${INDEL1} --known-sites ${INDEL2} --known-sites ${YBROWSE} -I ${SAMPLE%%.srt.bam}.bqsr.bam -O ${SAMPLE}.${CHR}.after -L ${CHR}
gatk-${GATK}/gatk AnalyzeCovariates --bqsr ${SAMPLE}.recal --before ${SAMPLE}.${CHR}.recal --after ${SAMPLE}.${CHR}.after --plots ${SAMPLE}.pdf
