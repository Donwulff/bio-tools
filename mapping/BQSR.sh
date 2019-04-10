#!/bin/sh

SAMPLE=${1:-sample1.srt.bam}
DATA=/mnt/GenomicData
CORES=`nproc`
GATK=4.1.1.0

# Script fetches and prepares necessary resources for BQSR and then drives a simple, single-machine BQSR workflow.
# Both GRCh37 and GRCh38 resources are fetched, although it might be better idea to liftover final results to GRCh37 if needed.

# Requires wget gawk tabix samtools bcftools python java

# BQSR is part of Broad Institute's Genome Analysis Toolkit Best Practices.
# It's a statistical method to empirically determine sequencing system error profile by masking out known variant sites.
# Because real novel variants are expected to distribute randomly on the reads, this gives a good idea of any systematic errors in sequencing.

# Explained: https://software.broadinstitute.org/gatk/documentation/article?id=11081
# Study: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048557/
# Known sites to use: https://software.broadinstitute.org/gatk/documentation/article.php?id=1247

# GATK Best Practices says to use "The most recent dbSNP release (build ID > 132)" for BaseRecalibrator.
# Some people prefer to stay with older snapshot, but I believe newer versions have much better coverage of rare populations.

# Standard resources don't cover Y chromosome well, so YBrowse SNP list is used for this.
# Best Practices excludes X and Y chromosomes, but you may wish to experiment what works best for your purposes.

if [ ! -e gatk-${GATK} ];
then
  wget -nc https://github.com/broadinstitute/gatk/releases/download/${GATK}/gatk-${GATK}.zip
  unzip gatk-${GATK}.zip
fi

set -x

# GRCh38 based with chr* chromosomes
REF=${DATA}/GRCh38/hg38.fa
DBSNP=${DATA}/GCF_000001405.38.dbSNP152.GRCh38p12.GATK.noINFO.vcf.gz
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

# Pre-prepared version of the dbSBP 152 GRCh38p12 file with all contigs that have them renamed to UCSC identifiers and INFO column blanked to minimize download and decompression
# This file is served off Sequencing.com
wget -nc https://api.sequencing.com/download.ashx?id=1c28f644-7c9d-43a1-aa68-2c2ee89a01da -O ${DATA}/GCF_000001405.38.dbSNP152.GRCh38p12.GATK.noINFO.vcf.gz
tabix ${DATA}/GCF_000001405.38.dbSNP152.GRCh38p12.GATK.noINFO.vcf.gz

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

  gawk -v RS="(\r)?\n" 'BEGIN { FS="\t" } !/^#/ { if ($10 != "na") print $7,$10; else print $7,$5 }' GCF_000001405.38_GRCh38.p12_assembly_report.txt > dbSNP-to-UCSC-GRCh38.p12.map
  time bcftools annotate --rename-chrs dbSNP-to-UCSC-GRCh38.p12.map GCF_000001405.38.bgz | gawk '/^#/ && !/^##contig=/ { print } !/^#/ { if( $1!="na" ) print }' \
    | bgzip -@${CORES} -l9 -c > ${DATA}/GCF_000001405.38.dbSNP152.GRCh38p12b.GATK.vcf.gz
  tabix ${DATA}/GCF_000001405.38.dbSNP152.GRCh38p12b.GATK.vcf.gz
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

# https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it
# "You can use the Mills_and_1000G_gold_standard.indels.hg38.vcf.gz and Homo_sapiens_assembly38.known_indels.vcf.gz as a replacement for the three original indel files."
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

# Assume dbSNP includes all contigs of interest
tabix -l ${DBSNP} | sort > dbSNP.contigs
samtools idxstats ${SAMPLE} | cut -f 1 | sort > ${SAMPLE}.contigs
join dbSNP.contigs ${SAMPLE}.contigs > recal.contigs

# TODO: Only run BaseRecalibrator on contigs we're going to use for the final calibration file; for now it can be used to check how they look.
if [ ! -e ${SAMPLE}.recal ];
then
  # Determine error profile: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php
  cat recal.contigs \
    | xargs -IZ -P${CORES} sh -c "
      if [ ! -e ${SAMPLE}.Z.recal ]; then \
      gatk-${GATK}/gatk --java-options -Xms4G BaseRecalibrator -R ${REF} \
      --known-sites ${DBSNP} \
      --known-sites ${INDEL1} \
      --known-sites ${INDEL2} \
      --known-sites ${YBROWSE} \
      -I ${SAMPLE} -O ${SAMPLE}.Z.recal -L Z ; \
      fi "

  # Tests show that the male X and Y chromosome covariates differ from the rest of the genome, perhaps due to being phased. (Test this hypothesis on female X later)
  # Now with YBrowse SNP's added, just use whole primary assembly covariates by default.
  ls ${SAMPLE}.chr[0-9].recal > chr.list
  ls ${SAMPLE}.chr[0-9][0-9].recal >> chr.list
  ls ${SAMPLE}.[0-9].recal >> chr.list
  ls ${SAMPLE}.[0-9][0-9].recal >> chr.list
  ls ${SAMPLE}.*X.recal ${SAMPLE}.*Y.recal >> chr.list
# Use all contigs with this instead
#  ls ${SAMPLE}.*.recal > chr.list
  gatk-${GATK}/gatk GatherBQSRReports -I chr.list -O ${SAMPLE}.recal

  # Generate a report showing difference between X and Y chromosome BQSR error covariates for checking.
  gatk-${GATK}/gatk AnalyzeCovariates --bqsr ${SAMPLE}.recal --before ${SAMPLE}.*X.recal --after ${SAMPLE}.*Y.recal --plots ${SAMPLE}.YvsX.pdf
fi

# Apply error profile: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php
# TODO: We can do this in paraller, although copying piecewise BAM's back together will require extra IO.
if [ ! -e ${SAMPLE%%.srt.bam}.bqsr.bam ];
then
  gatk-${GATK}/gatk ApplyBQSR -R ${REF} --bqsr ${SAMPLE}.recal -I ${SAMPLE} -O ${SAMPLE%%.srt.bam}.bqsr.bam
fi

# Check for residual error: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_AnalyzeCovariates.php
CHR=`tabix -l ${DBSNP} | grep -m1 '20$'`
gatk-${GATK}/gatk --java-options -Xms4G BaseRecalibrator -R ${REF} \
  --known-sites ${DBSNP} --known-sites ${INDEL1} --known-sites ${INDEL2} --known-sites ${YBROWSE} -I ${SAMPLE%%.srt.bam}.bqsr.bam -O ${SAMPLE}.${CHR}.after -L ${CHR}
gatk-${GATK}/gatk AnalyzeCovariates --bqsr ${SAMPLE}.recal --before ${SAMPLE}.${CHR}.recal --after ${SAMPLE}.${CHR}.after --plots ${SAMPLE}.pdf
