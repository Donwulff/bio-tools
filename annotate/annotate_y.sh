#!/bin/sh

# This script was used to annotate Y-features on the Tyrolean Iceman sequence
# At the time GRCh38 mapped features weren't yet available, so it does liftover
# Saving this one for reference, it shows how to liftover and annotate when needed

# Build programs
cd src
wget http://hgdownload.cse.ucsc.edu/admin/exe/userApps.v334.src.tgz
tar -zxvvf userApps.v334.src.tgz
cd userApps
make

# These aren't actually used right now, they were out of date so I dumped features straight out of gbrowse
wget http://ybrowse.org/gbrowse2/gff/snps_hg19.gff3
wget http://ybrowse.org/gbrowse2/gff/snps_hg38.gff3

# Change capped Chr to chr. STR likely can't be annotated, but do them just as well
cat /mnt/gbrowse2/snps.gff3 /mnt/gbrowse2/indel.gff3 /mnt/gbrowse2/STR.gff3 | sed s/^Chr/chr/g > snps.gff3

# Get the chain file for this conversion
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Run the conversion
bin/liftOver -gff snps.gff3 hg19ToHg38.over.chain.gz snpsHg38.gff3 snpsHg38.log

# It tells you to run this instead, but it generates empty gp file!
bin/ldHgGene -out=snpsHg37.gp phylo snps snps.gff3

# This might magically sort #'s to start, then number order for tabix file
sort -k4n snpsHg38.gff3 | bgzip > snpsHg38.gff3.gz
tabix snpsHg38.gff3.gz

# Header line to add to vcf needed to lift up the annotation
echo '##INFO=<ID=Name,Number=1,Type=String,Description="Yet another header line">' > header.add

# bcftools annotate is broken, can't lift INFO to INFO, so just sed the broken Name part out.
bcftools annotate -a snpsHg38.gff.gz -h header.add -c CHROM,-,-,FROM,TO,-,-,-,INFO/Name /mnt/Iceman/GRC13292546_bowtie2_sorted_dedupe_numeric_annotate.vcf.gz \
  | sed s/Name=Name/Name/g | bgzip > /mnt/Iceman/GRC13292546_bowtie2_sorted_dedupe_numeric_annotate_phylog.vcf.gz
