#!/bin/sh

# This simple script constructs an EXPERIMENTAL bwa alignment index including hs38s1 decoy sequences, released patches and HLA sequences.
# The genomic coordinates will stay same with GRCh38, but many reads will add preferentially to the added contigs.
# There's probably no good reason one should use this reference - it isn't compatible with anything, and continues to change.
# Bioinformatics relies on stable, comparable results between runs, and there aren't existing tools that would use the additional contigs.
# However, if you're willing to experiment and don't care about results that are statistically comparable with other versions, go ahead.

# Studies have shown that the inclusion of the decoy sequences improves mapping results, and the alternate sequences conceivably do the same.
# See explanation & preliminary results: https://github.com/lh3/bwa/blob/master/README-alt.md

# Some more general references:
# http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use explains the choice of the base reference.
# Original hs37d5 decoy sequence presentation http://lh3lh3.users.sourceforge.net/download/decoyseq.pdf
# Some notes in ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707
# Scripts mentioned are available in https://github.com/lh3/misc/tree/master/seq/novoseq

# hs38d1 construction is covered in https://www.ncbi.nlm.nih.gov/pubmed/27654912 Supplementary section 5
# Worked example of alt contigs: https://software.broadinstitute.org/gatk/documentation/article?id=8017

# Naming is a complicated issue. Currently I'm using UCSC naming only for hg38Patch11, with patch 12 and 13 in NCBI because they're not official yet.
# Patch 12 whole genome is now available in separate folder, but isn't masked for analysis use, and isn't available in patches-only format.
# HLA contigs have been renamed to use the HLA allele number notation, because many tools can't handle special characters, and the numbers don't change.

# Are any of the patches for regions that are masked in the analysis set? Alt contigs shouldn't break things, but those would be spurious.

set -x

# WARNING! Oral microbiome is experimental.
VERSION_BASE="hg38"
VERSION_PATCH="p13"
VERSION_DECOY="D"
VERSION_HLA="-"
VERSION_ORAL="O911"
VERSION_EXTRA="alt"
# VERSION_HLA will be read from IMGT/HLA Allele_status.txt file upon import of hla_gen.fasta

# National Center for Biotechnology Information Analysis Set https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#seqsforalign
# We currently need hs381d1 with UCSC naming, and the assembly with PAR & centromeric masking. We could get these from full_plus_hs38d1 BWA index, but lets download.
[ -e GCA_000786075.2_hs38d1_genomic_unmapped.alt ] || \
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

# BWA alignment set without alt contigs or decoy sequences is used to determine mapping of alt contigs into the primary assembly in alt file
if [ ! -e additional_hg38_contigs.alt ] && [ ! -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sa ];
then
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
  tar zkxf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
fi

# We need this only for the base ALT-file, perhaps we could regenerate it ourselves; there's also copy in the git repo to skip expensive calculations
if [ ! -f GCA_000001405.15_GRCh38_full_analysis_set.fna.alt ];
then
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz
  tar zkxf GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz
fi

# University of California Santa Cruz UCSC's contig names used in their Golden Path genome browser have become standard http://genome.ucsc.edu/
wget -nc http://hgdownload.cse.ucsc.edu/goldenPath/hg38/hg38Patch11/hg38Patch11.fa.gz
# p12 UCSC names aren't yet in the p12 genome assembly report, so the names in this version are used in future; ignore currently
#wget -nc ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz

# Get the NCBI reference genomes so we can construct incremental patches between them to add to the UCSC reference genome
if [ ! -e GRCh38Patch12.fa.gz ] || [ ! -e GRCh38Patch13.fa.gz ]; then
  # Genome Reference Consortium https://www.ncbi.nlm.nih.gov/grc/human releases cumulative patches to the latest assembly
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.26_GRCh38.p11/GCA_000001405.26_GRCh38.p11_genomic.fna.gz
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GCA_000001405.27_GRCh38.p12_genomic.fna.gz
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz

  # I couldn't find a source for incremental patches to the human assembly, so we need to diff and clean it up.
  # These could be converted to the UCSC naming, but because they're not yet officially in UCSC, that could be misleading.
  # This is done stepwise to keep the new patches at the end of the file in case one needs to compare mappings etc.
  # Patch 13 includes some changes that do not effect reference, ignore case and contig state change.
  zdiff --ignore-case GCA_000001405.26_GRCh38.p11_genomic.fna.gz GCA_000001405.27_GRCh38.p12_genomic.fna.gz | grep "^> " | grep -v "UNVERIFIED_ORG" | cut -c3- | gzip -c > GRCh38Patch12.fa.gz
  zdiff --ignore-case GCA_000001405.27_GRCh38.p12_genomic.fna.gz GCA_000001405.28_GRCh38.p13_genomic.fna.gz | grep "^> " | grep -v "UNVERIFIED_ORG" | cut -c3- | gzip -c > GRCh38Patch13.fa.gz
fi

# European Molecular Biology Laboratory publishes the IPD-IMGT/HLA database with World Health Organization's naming https://www.ebi.ac.uk/ipd/imgt/hla/ nb. this DOES change a lot
# To regenerate, delete Allele_status.txt hla_gen.fasta
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allele_status.txt
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta

VERSION_HLA="H$(grep version Allele_status.txt | tr -cd '[0-9]')"
VERSION="$VERSION_BASE$VERSION_PATCH$VERSION_DECOY$VERSION_HLA$VERSION_ORAL$VERSION_EXTRA"

# Convert the HLA FASTA sequence names and compress it, no longer using bwa-kit HLA allele notation because : and * mess up most tools!
[ -e hla_gen.$VERSION_HLA.fasta.gz ] || sed "s/^>HLA:/>/" hla_gen.fasta | gzip -c > hla_gen.$VERSION_HLA.fasta.gz

# Construct mapping index for whole assembly + HLA to compare decoys and microbiome against
if [ ! -e hg38.p12.p13.$VERSION_HLA.fa.gz.sa ]; then
  cat GCA_000001405.15_GRCh38_full_analysis_set.fna.gz \
      hg38Patch11.fa.gz \
      GRCh38Patch12.fa.gz \
      GRCh38Patch13.fa.gz \
      hla_gen.$VERSION_HLA.fasta.gz \
    > hg38.p12.p13.$VERSION_HLA.fa.gz
  bwa index hg38.p12.p13.$VERSION_HLA.fa.gz
fi

## Steps to clean up decoy sequences
DECOY_BASE=GCA_000786075.2_hs38d1_p13_${VERSION_HLA}_genomic
if [ "$VERSION_DECOY" != "" ] && [ ! -e ${DECOY_BASE}_unmapped.alt ]; then
  # Filter out decoys which map to the current assembly for 101bp or more
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz
  bwa mem -t`nproc` -k101 hg38.p12.p13.$VERSION_HLA.fa.gz GCA_000786075.2_hs38d1_genomic.fna.gz > $DECOY_BASE.sam

  # Rename unmapped decoy contigs into the UCSC style used by reference genomes
  samtools view -f0x4 $DECOY_BASE.sam | \
    gawk -v OFS="\t" '{ gsub("\\.","v",$1); print "chrUn_"$1"_decoy"; }' > \
      ${DECOY_BASE}_unmapped.list

  # Use the unmapped contigs to select matching decoy sequences from the analysis set to easily rename them and remove soft-masking
  [ -e GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ] || gzip -cd GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz > \
    GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
  samtools faidx GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
  samtools faidx GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna -r ${DECOY_BASE}_unmapped.list | \
    bgzip -c > ${DECOY_BASE}_unmapped.fna.gz

  # Generate alt lines from alignment of remaining, unmapped decoys against the primary assy
  samtools view -f0x4 $DECOY_BASE.sam | \
    gawk -v OFS="\t" '{ gsub("\\.","v",$1); $1 = "chrUn_"$1"_decoy"; $6 = length($10)"M"; $10 = "*"; $11 = "*"; NF=11; print }' > \
      ${DECOY_BASE}_unmapped.alt
fi

## The Forsyth "expanded Human Oral Microbiome Database" http://www.homd.org
ORAL_BASE=oral_microbiome_p13_${VERSION_HLA}_${VERSION_ORAL}
if [ "$VERSION_ORAL" != "" ] && [ ! -e ${ORAL_BASE}_unmapped.alt ]; then
  # Filter out decoys which map to the current assembly for 101bp or more
  wget -nc http://www.homd.org/ftp/genomes/PROKKA/V9.11/fsa/ALL_genomes.fsa
  bwa mem -t`nproc` -k101 hg38.p12.p13.$VERSION_HLA.fa.gz ALL_genomes.fsa > $ORAL_BASE.sam
  samtools view -f0x4 $ORAL_BASE.sam | cut -f1 > \
    ${ORAL_BASE}_unmapped.list

  # Use the unmapped contigs to select matching decoy sequences from the analysis set
  samtools faidx ALL_genomes.fsa
  samtools faidx ALL_genomes.fsa -r ${ORAL_BASE}_unmapped.list | \
    bgzip -c > ${ORAL_BASE}_unmapped.fna.gz

  # Generate alt lines from alignment of remaining, unmapped decoys against the primary assy
  samtools view -f0x4 $ORAL_BASE.sam | \
    gawk -v OFS="\t" '{ $6 = length($10)"M"; $10 = "*"; $11 = "*"; NF=11; print }' > \
      ${ORAL_BASE}_unmapped.alt
fi

# Scoring parameters found counting Alignment Score from bwakit hg38DH.fa.alt; this generates more supplementary alignments and missed odd MapQ 30 line
if [ ! -e additional_hg38_contigs.alt ]; then
  cat hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz hla_gen.$VERSION_HLA.fasta.gz > additional_hg38_p13_${VERSION_HLA}_contigs.fa.gz
  bwa mem -t`nproc` -A2 -B3 -O4 -E1 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna additional_hg38_p13_${VERSION_HLA}_contigs.fa.gz \
    | samtools view -q60 - \
    | gawk '{ OFS="\t"; $10 = "*"; print }' > additional_hg38_p13_${VERSION_HLA}_contigs.alt
fi

zcat GCA_000001405.15_GRCh38_full_analysis_set.fna.gz ${DECOY_BASE}_unmapped.fna.gz \
     hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz hla_gen.$VERSION_HLA.fasta.gz ${ORAL_BASE}_unmapped.fna.gz > ${VERSION}.fa
cat GCA_000001405.15_GRCh38_full_analysis_set.fna.alt \
    ${DECOY_BASE}_unmapped.alt \
    additional_hg38_p13_${VERSION_HLA}_contigs.alt \
    ${ORAL_BASE}_unmapped.alt \
  > ${VERSION}.fa.alt

bwa index ${VERSION}.fa

samtools faidx ${VERSION}.fa
samtools dict -a "GRCh38" -s "Homo Sapiens" -u "${VERSION}.fa" ${VERSION}.fa -o ${VERSION}.dict

#gatk-4.1.2.0/gatk FindBadGenomicKmersSpark -R ${VERSION}.fa -O ${VERSION}.fa.txt
