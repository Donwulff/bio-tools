#!/bin/sh

# This simple script constructs an EXPERIMENTAL bwa alignment index including hs38s1 decoy sequences, released patches and HLA sequences.
# There's probably no good reason one should use this reference - it isn't compatible with anything, and continues to change.
# Bioinformatics relies on stable, comparable results between runs, and there aren't existing tools that would use the additional contigs.
# Studies have shown that the inclusion of the decoy sequences improves mapping results, and the alternate sequences might do the same.

# http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use explains the choice of the base reference.
# Decoy sequence presentation http://lh3lh3.users.sourceforge.net/download/decoyseq.pdf

# TODO: Repeat parts of the decoy sequence construction to ensure there's no interaction with the added contigs
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707
# After that it may be possible to regenerate the decoy sequences from desired sources and versions

# Some of the HLA should overlap with the 35 HLA versions in the reference. Do we need to remove duplicates?
# Are any of the patches for regions masked in the analysis set? Alt contigs shouldn't break things, but they're spurious.

# This contains the bwa alt contig definition file, but we recreate the indexes, so they're not needed.
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz
#wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/unmasked_cognates_of_masked_CEN_PAR.txt

# National Center for Biotechnology Information Analysis Set https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#seqsforalign
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz

# University of California Santa Cruz UCSC's contig names used in their Golden Path genome browser have become standard http://genome.ucsc.edu/
wget -nc http://hgdownload.cse.ucsc.edu/goldenPath/hg38/hg38Patch11/hg38Patch11.fa.gz

# Genome Reference Consortium https://www.ncbi.nlm.nih.gov/grc/human releases cumulative patches to the latest assembly
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GCA_000001405.26_GRCh38.p11_genomic.fna.gz
wget -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GCA_000001405.27_GRCh38.p12_genomic.fna.gz

# European Molecular Biology Laboratory publishes the IPD-IMGT/HLA database with World Health Organization's naming https://www.ebi.ac.uk/ipd/imgt/hla/
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta

# I couldn't find a source for incremental patches to the human assembly, so we need to diff and clean it up.
# These could be converted to the UCSC naming, but because they're not yet officially in UCSC, that could be misleading.
zdiff GCA_000001405.26_GRCh38.p11_genomic.fna.gz GCA_000001405.27_GRCh38.p12_genomic.fna.gz | grep "^> " | cut -c3- | gzip -c > GRCh38patch12.fa.gz

# Convert the HLA FASTA sequence names into the format used by bwa's bwa-kit release and compress it
sed "s/^>HLA:HLA..... />HLA-/" hla_gen.fasta | gzip -c > hla_gen.fasta.gz

# There's no documentation for how the bwa alt-file should be constructed. This is just a basic starting point.
# https://github.com/lh3/bwa/blob/master/README-alt.md
bwa mem -x intractg -t4 hs38.fa hg38Patch11.fa.gz \
  | samtools view - \
  | gawk '{ OFS="\t"; $10 = "*"; print }' > hg38Patch11.map

bwa mem -x intractg -t4 hs38.fa GRCh38patch12.fa.gz \
  | samtools view - \
  | gawk '{ OFS="\t"; $10 = "*"; print }' > GRCh38patch12.map

bwa mem -x intractg -t4 hs38.fa hla_gen.fasta.gz \
  | samtools view - \
  | gawk '{ OFS="\t"; $10 = "*"; print }' > hla_gen.map

zcat GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz hg38Patch11.fa.gz GRCh38patch12.fa.gz hla_gen.fasta.gz > hg38-all.fa
cat GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.alt hg38Patch11.map GRCh38patch12.map hla_gen.map > hg38-all.fa.alt

bwa index hg38-all.fa
