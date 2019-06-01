#!/bin/sh

# This simple script constructs an EXPERIMENTAL bwa alignment index including hs38s1 decoy sequences, released patches and HLA sequences.
# There's probably no good reason one should use this reference - it isn't compatible with anything, and continues to change.
# Bioinformatics relies on stable, comparable results between runs, and there aren't existing tools that would use the additional contigs.
# Studies have shown that the inclusion of the decoy sequences improves mapping results, and the alternate sequences might do the same.

# http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use explains the choice of the base reference.
# Decoy sequence presentation http://lh3lh3.users.sourceforge.net/download/decoyseq.pdf
# hs38d1 is actually detailed in https://www.ncbi.nlm.nih.gov/pubmed/27654912 Supplementary section 5

# Are any of the patches for regions masked in the analysis set? Alt contigs shouldn't break things, but they're spurious.

# BWA alignment set without alt contigs or decoy sequences is used to determine mapping of alt contigs into the primary assembly for alt file
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
if [ ! -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sa ];
then
  tar zkxf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
fi

# National Center for Biotechnology Information Analysis Set https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#seqsforalign
# We currently need hs381d1 with UCSC naming, and the assembly with PAR & centromeric matching. Possible to get these otherwise?
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.fai
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

# University of California Santa Cruz UCSC's contig names used in their Golden Path genome browser have become standard http://genome.ucsc.edu/
wget -nc http://hgdownload.cse.ucsc.edu/goldenPath/hg38/hg38Patch11/hg38Patch11.fa.gz
# p12 UCSC names aren't yet in the p12 genome assembly report, so the names in this version are future
wget -nc ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz

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
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta

# Convert the HLA FASTA sequence names into the format used by bwa's bwa-kit release and compress it
sed "s/^>HLA:HLA..... />HLA-/" hla_gen.fasta | gzip -c > hla_gen.fasta.gz

## Steps to clean up decoy sequences
if [ -e GCA_000786075.2_hs38d1_genomic_unmapped.fna.gz ]; then
  # Construct mapping index for whole assembly + HLA
  if [ -e hg38.p12.p13.hla.fa.gz ]; then
    cat hg38.p12.fa.gz GRCh38Patch13.fa.gz hla_gen.fasta.gz > hg38.p12.p13.hla.fa.gz
    bwa index hg38.p12.p13.hla.fa.gz
  fi

  # Filter out decoys which map to the current assembly, and re-generate alt lines
  wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz
  bwa mem -k101 hg38.p12.p13.hla.fa.gz GCA_000786075.2_hs38d1_genomic.fna.gz > GCA_000786075.2_hs38d1_genomic.sam
  samtools view -f0x4 GCA_000786075.2_hs38d1_genomic.sam | \
    gawk -v OFS="\t" '{ gsub("\\.","v",$1); $1 = "chrUn_"$1"_decoy"; $6 = length($10)"M"; $10 = "*"; $11 = "*"; NF=11; print }' \
    GCA_000786075.2_hs38d1_genomic_unmapped.sam

  # Use the alt-lines to select matching decoy sequences from the analysis set
  cut -f1 GCA_000786075.2_hs38d1_genomic_unmapped.sam > GCA_000786075.2_hs38d1_genomic.list
  samtools faidx GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
  samtools faidx GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz -r GCA_000786075.2_hs38d1_genomic.list | \
    bgzip -c > GCA_000786075.2_hs38d1_genomic_unmapped.fna.gz
fi

# There's no documentation for how the bwa alt-file should be constructed. This is just a basic starting point.
# https://github.com/lh3/bwa/blob/master/README-alt.md
cat hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz hla_gen.fasta.gz > additional_hg38_contigs.fa.gz
bwa mem -t`nproc`-x intractg -t4 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna additional_hg38_contigs.fa.gz \
  | samtools view - \
  | gawk '{ OFS="\t"; $10 = "*"; print }' > additional_hg38_contigs.map

zcat GCA_000001405.15_GRCh38_full_analysis_set.fna.gz GCA_000786075.2_hs38d1_genomic_unmapped.fna.gz \
     hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz hla_gen.fasta.gz > hg38DH-p13.fa
cat GCA_000001405.15_GRCh38_full_analysis_set.fna.alt GCA_000786075.2_hs38d1_genomic_unmapped.sam additional_hg38_contigs.map > hg38DH-p13.fa.alt

bwa index hg38DH-p13.fa

#gatk-4.1.2.0/gatk FindBadGenomicKmersSpark -R hg38DH-p13.fa -O hg38DH-p13.fa.txt
