#!/bin/sh

# This simple script constructs an EXPERIMENTAL bwa alignment index including hs38s1 decoy sequences, released patches and HLA sequences.
# The genomic coordinates will stay same with GRCh38, but many reads will add preferentially to the added contigs.
# There's probably no good reason one should use this reference - it isn't compatible with anything, and continues to change.
# Bioinformatics relies on stable, comparable results between runs, and there aren't existing tools that would use the additional contigs.
# However, if you're willing to experiment and don't care about results that are statistically comparable with other versions, go ahead.

# Studies have shown that the inclusion of the decoy sequences improves mapping results, and the alternate sequences conceivably do the same.
# See explanation & preliminary results: https://github.com/lh3/bwa/blob/master/README-alt.md

# Some more general references:
# https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use explains the choice of the base reference.
# Original hs37d5 decoy sequence presentation https://lh3lh3.users.sourceforge.net/download/decoyseq.pdf
# Some notes in https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707
# Scripts mentioned are available in https://github.com/lh3/misc/tree/master/seq/novoseq

# hs38d1 construction is covered in https://www.ncbi.nlm.nih.gov/pubmed/27654912 Supplementary section 5
# Worked example of alt contigs: https://software.broadinstitute.org/gatk/documentation/article?id=8017

# Naming is a complicated issue. Currently I'm using UCSC naming only for hg38Patch11, with patch 12 and 13 in NCBI because they're not official yet.
# Patch 12 whole genome is now available in separate folder, but isn't masked for analysis use, and isn't available in patches-only format.
# HLA contigs have been renamed to use the HLA allele number notation, because many tools can't handle special characters, and the numbers don't change.

# Are any of the patches for regions that are masked in the analysis set? Alt contigs shouldn't break things, but those would be spurious.

set -x

# Changing VERSION_PATCH doesn't affect built reference in itself, it's just for reference.
# VERSION_HLA chooses between HLA 'H' and accession number 'N' nomenclature.
# HLA alleles contain asterisks and other special letters which may impede downstream processing.
# VERSION_HLA will be read from latest IPD-IMGT/HLA data file,
# VERSION_DECOY '' means not to include decoy sequences,
# VERSION_EXTRA 'hla' means to use HLA-* allele names which cause issues with some downstream analysis.
# WARNING! Oral microbiome is experimental. Version 10.01 is huge, and requites 64GB memory; 9.15 will fit in 32GB memory.
VERSION_BASE="hg38"
VERSION_PATCH="p14"
VERSION_DECOY="D"
VERSION_HLA="A"
VERSION_ORAL="10.1"
VERSION_EXTRA=""

# Extra eXperimental build chm13 reference instead.
#VERSION_BASE="chm13"
#VERSION_PATCH="v2.0_maskedY_rCRS"

# National Center for Biotechnology Information Analysis Set https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#seqsforalign
# We currently need hs381d1 with UCSC naming, the assembly with PAR & centromeric masking, and no-alt set for masked alt-alignment creation.
# It's starting to feel like we would be better off creating these from full_plus_hs38d1 BWA index, but lets download.
if [ "$VERSION_BASE" = "hg38" ]; then
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
fi

# We need this only for the base ALT-file, perhaps we could regenerate it ourselves; there's also copy in the git repo to skip expensive calculations
if [ "$VERSION_BASE" = "hg38" ] && [ ! -f GCA_000001405.15_GRCh38_full_analysis_set.fna.alt ];
then
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz
  tar zkxf GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz
fi

# Validate downloaded files; we should do this before extracting, but for simplicity, fix these manually if needed.
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/md5checksums.txt
rm GCA_000001405.15_GRCh38_GRC_exclusions.bed
md5sum --ignore-missing --check md5checksums.txt
if [ "$VERSION_BASE" = "hg38" ] && [ "$?" != "0" ]; then
  echo "Problem with downloading base references."
  exit 1
fi

# Funny thing, the md5checksums.txt has not been updated after p14 change of the GCA_000001405.15_GRCh38_GRC_exclusions.bed, so we can't check it.
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed

# University of California Santa Cruz UCSC's contig names used in their Golden Path genome browser have become standard https://genome.ucsc.edu/
[ "$VERSION_BASE" = "hg38" ] && wget -nc https://hgdownload.cse.ucsc.edu/goldenPath/hg38/hg38Patch11/hg38Patch11.fa.gz
# p12 UCSC names aren't yet in the p12 genome assembly report, so the names in this version are used in future; ignore currently
#wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz

# Get the NCBI reference genomes so we can construct incremental patches between them to add to the UCSC reference genome
if [ "$VERSION_BASE" = "hg38" ] && ( [ ! -e GRCh38Patch12.fa.gz ] || [ ! -e GRCh38Patch13.fa.gz ] || [ ! -e GRCh38Patch14.fa.gz ] ); then
  # Genome Reference Consortium https://www.ncbi.nlm.nih.gov/grc/human releases cumulative patches to the latest assembly
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.26_GRCh38.p11/GCA_000001405.26_GRCh38.p11_genomic.fna.gz
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GCA_000001405.27_GRCh38.p12_genomic.fna.gz
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.26_GRCh38.p11/md5checksums.txt -O md5checksums.p11.txt
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/md5checksums.txt -O md5checksums.p12.txt
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/md5checksums.txt -O md5checksums.p13.txt
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/md5checksums.txt -O md5checksums.p14.txt
  cat md5checksums.p11.txt md5checksums.p12.txt md5checksums.p13.txt md5checksums.p14.txt > md5checksums.all.txt
  md5sum --ignore-missing --check md5checksums.all.txt
  if [ "$?" != "0" ]; then
    echo "Problem with downloading patch references for differential."
    exit 1
  fi

  # I couldn't find a source for incremental patches to the human assembly, so we need to diff and clean it up.
  # These could be converted to the UCSC naming, but because they're not yet officially in UCSC, that could be misleading.
  # This is done stepwise to keep the new patches at the end of the file in case one needs to compare mappings etc.
  # Patch 13 includes some changes that do not affect sequence, ignore case and contig state change.
  samtools dict -a GRCh38 -s "Homo Sapiens" -u "" GCA_000001405.26_GRCh38.p11_genomic.fna.gz -o GCA_000001405.26_GRCh38.p11_genomic.fna.dict
  samtools dict -a GRCh38 -s "Homo Sapiens" -u "" GCA_000001405.27_GRCh38.p12_genomic.fna.gz -o GCA_000001405.27_GRCh38.p12_genomic.fna.dict
  samtools dict -a GRCh38 -s "Homo Sapiens" -u "" GCA_000001405.28_GRCh38.p13_genomic.fna.gz -o GCA_000001405.28_GRCh38.p13_genomic.fna.dict
  samtools dict -a GRCh38 -s "Homo Sapiens" -u "" GCA_000001405.29_GRCh38.p14_genomic.fna.gz -o GCA_000001405.29_GRCh38.p14_genomic.fna.dict

  gzip -cd GCA_000001405.29_GRCh38.p14_genomic.fna.gz | bgzip > GCA_000001405.29_GRCh38.p14_genomic.bgzip.fna.gz
  samtools faidx GCA_000001405.29_GRCh38.p14_genomic.bgzip.fna.gz

  # GRCh38Patch14 is different, KQ983257.1 and KQ983258.1 are in different order from previous patches, whhich makes KQ983257.1 look like a new sequence; keeping original ordering.
  # KQ759759.1 and KQ759762.1 from Patch12 are replaced with new, revision 2 sequences. They could be left out, but for downstream analysis I'm keeping them.
  # On 20th April 2022, location of the contig KI270825.1 was moved to different place in the .p14 reference, in which it was already declared contamination.
  diff -u0 GCA_000001405.26_GRCh38.p11_genomic.fna.dict GCA_000001405.27_GRCh38.p12_genomic.fna.dict | grep "^+@SQ" | cut -f2 | cut -d':' -f2 > GRCh38Patch12.list
  diff -u0 GCA_000001405.27_GRCh38.p12_genomic.fna.dict GCA_000001405.28_GRCh38.p13_genomic.fna.dict | grep "^+@SQ" | cut -f2 | cut -d':' -f2 > GRCh38Patch13.list
  diff -u0 GCA_000001405.28_GRCh38.p13_genomic.fna.dict GCA_000001405.29_GRCh38.p14_genomic.fna.dict | grep "^+@SQ" | cut -f2 | cut -d':' -f2 | grep -v "KQ983257\.1" | grep -v "KI270825\.1" > GRCh38Patch14.list

  samtools faidx GCA_000001405.29_GRCh38.p14_genomic.bgzip.fna.gz -r GRCh38Patch12.list | gzip -c > GRCh38Patch12.fa.gz
  samtools faidx GCA_000001405.29_GRCh38.p14_genomic.bgzip.fna.gz -r GRCh38Patch13.list | gzip -c > GRCh38Patch13.fa.gz
  samtools faidx GCA_000001405.29_GRCh38.p14_genomic.bgzip.fna.gz -r GRCh38Patch14.list | gzip -c > GRCh38Patch14.fa.gz
fi

# European Molecular Biology Laboratory publishes the IPD-IMGT/HLA database with World Health Organization's naming https://www.ebi.ac.uk/ipd/imgt/hla/ nb. this DOES change a lot
# To regenerate with latest (possibly updated) version, delete Allele_status.txt hla_gen.fasta
wget -nc https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allele_status.txt
wget -nc https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta

# Convert the HLA FASTA sequence names and compress it, no longer using bwa-kit HLA allele notation by default because : and * mess up most tools!
if [ "$VERSION_HLA" = "A" ]; then
  VERSION_HLA="${VERSION_HLA}$(grep version Allele_status.txt | tr -cd '[0-9]')"
  [ -e "hla_gen.${VERSION_HLA}.fasta.gz" ] || sed "s/^>HLA:HLA[0-9]* />HLA-/" hla_gen.fasta | gzip -c > hla_gen.${VERSION_HLA}.fasta.gz
else
  VERSION_HLA="$VERSION_HLA$(grep version Allele_status.txt | tr -cd '[0-9]')"
  [ -e "hla_gen.${VERSION_HLA}.fasta.gz" ] || sed "s/^>HLA:/>/" hla_gen.fasta | gzip -c > hla_gen.${VERSION_HLA}.fasta.gz
fi

VERSION_ORAL_CODE=${VERSION_ORAL_URL//./}
VERSION=${VERSION_BASE}${VERSION_PATCH}${VERSION_DECOY}${VERSION_HLA}${VERSION_ORAL_CODE}${VERSION_EXTRA}

# Construct mapping index for whole assembly + HLA to compare decoys and microbiome against
if [ "${VERSION_BASE}" = "hg38" ] && [ ! -e "${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz.sa" ]; then
  if [ ! -e "GCA_000001405.15_GRCh38_full_analysis_set_masked.fna.gz" ]; then
    gzip -cd GCA_000001405.15_GRCh38_full_analysis_set.fna.gz > GCA_000001405.15_GRCh38_full_analysis_set.fna
    bedtools maskfasta -fullHeader -fi GCA_000001405.15_GRCh38_full_analysis_set.fna -fo GCA_000001405.15_GRCh38_full_analysis_set_masked.fna -bed GCA_000001405.15_GRCh38_GRC_exclusions.bed
    gzip GCA_000001405.15_GRCh38_full_analysis_set_masked.fna
  fi
  # Use concatenated gzip's for speed because this is temporary index.
  cat GCA_000001405.15_GRCh38_full_analysis_set_masked.fna.gz \
      hg38Patch11.fa.gz \
      GRCh38Patch12.fa.gz \
      GRCh38Patch13.fa.gz \
      GRCh38Patch14.fa.gz \
      hla_gen.${VERSION_HLA}.fasta.gz \
    > ${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz
  bwa index ${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz
fi

if [ "${VERSION_BASE}" = "chm13" ] && [ ! -e "${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz.sa" ]; then
  wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/${VERSION_BASE}${VERSION_PATCH}.fa.gz
  # Use concatenated gzip's for speed because this is temporary index.
  cat ${VERSION_BASE}${VERSION_PATCH}.fa.gz \
      hla_gen.${VERSION_HLA}.fasta.gz \
    > ${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz
  bwa index ${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz
fi

## Steps to clean up decoy sequences
DECOY_BASE=GCA_000786075.2_hs38d1_${VERSION_BASE}${VERSION_PATCH}_${VERSION_HLA}_genomic
if [ "${VERSION_DECOY}" != "" ] && [ ! -e "${DECOY_BASE}_unmapped.alt" ]; then
  # Filter out decoys which map to the current assembly for 101bp or more
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz
  bwa mem -t`nproc` -k101 ${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz GCA_000786075.2_hs38d1_genomic.fna.gz > ${DECOY_BASE}.sam

  if [ "$VERSION_BASE" = "hg38" ]; then
    # Rename unmapped decoy contigs into the UCSC style used by reference genomes
    samtools view -f0x4 ${DECOY_BASE}.sam | \
      gawk -v OFS="\t" '{ gsub("\\.","v",$1); print "chrUn_"$1"_decoy"; }' > \
        ${DECOY_BASE}_unmapped.list

    # Use the unmapped contigs to select matching decoy sequences from the analysis set to easily rename them and remove soft-masking
    [ -e GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ] || gzip -cd GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz > \
      GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
    samtools faidx GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
    samtools faidx GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna -r ${DECOY_BASE}_unmapped.list | \
      bgzip -c > ${DECOY_BASE}_unmapped.fna.gz
  else
    # In chm13 we don't need the Analysis Ready reference genome, so let's extract the decoy sequences with original names
    samtools view -f0x4 ${DECOY_BASE}.sam | \
      cut -f1 > \
      ${DECOY_BASE}_unmapped.list

    # This preservers lowe-case soft-masking; bwa itself won't use it, though we really should use the chm13 masks for it otherwise
    [ -e GCA_000786075.2_hs38d1_genomic.fna ] || gzip -cd GCA_000786075.2_hs38d1_genomic.fna.gz > \
      GCA_000786075.2_hs38d1_genomic.fna
    samtools faidx GCA_000786075.2_hs38d1_genomic.fna
    samtools faidx GCA_000786075.2_hs38d1_genomic.fna -r ${DECOY_BASE}_unmapped.list | \
      gawk -v OFS="\t" '/^>/ { gsub("\\.","v",$1); gsub("^>","",$1); $1=">chrUn_"$1"_decoy"; } { print; }' | \
      bgzip -c > ${DECOY_BASE}_unmapped.fna.gz
  fi

  # Generate alt lines from alignment of remaining, unmapped decoys against the primary assy
  samtools view -f0x4 ${DECOY_BASE}.sam | \
    gawk -v OFS="\t" '{ gsub("\\.","v",$1); $1 = "chrUn_"$1"_decoy"; $6 = length($10)"M"; $10 = "*"; $11 = "*"; NF=11; print }' > \
      ${DECOY_BASE}_unmapped.alt
fi

## The Forsyth "expanded Human Oral Microbiome Database" https://www.homd.org
ORAL_BASE=oral_microbiome_${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}_${VERSION_ORAL_CODE}_genomic
if [ "${VERSION_ORAL}" != "" ] && [ ! -e "${ORAL_BASE}_unmapped.alt" ]; then
  # Filter out decoys which map to the current assembly for 101bp or more
  wget -nc https://www.homd.org/ftp/genomes/PROKKA/V${VERSION_ORAL}/fsa/ALL_genomes.fsa -O ${ORAL_BASE}.fsa
  bwa mem -t`nproc` -k101 ${VERSION_BASE}${VERSION_PATCH}.${VERSION_HLA}.fa.gz ${ORAL_BASE}.fsa > ${ORAL_BASE}.sam
  samtools view -f0x4 ${ORAL_BASE}.sam | cut -f1 > \
    ${ORAL_BASE}_unmapped.list

  # Use the unmapped contigs to select matching decoy sequences from the analysis set
  samtools faidx ${ORAL_BASE}.fsa
  samtools faidx ${ORAL_BASE}.fsa -r ${ORAL_BASE}_unmapped.list | \
    bgzip -c > ${ORAL_BASE}_unmapped.fna.gz

  # Generate alt lines from alignment of remaining, unmapped decoys against the primary assy
  samtools view -f0x4 ${ORAL_BASE}.sam | \
    gawk -v OFS="\t" '{ $6 = length($10)"M"; $10 = "*"; $11 = "*"; NF=11; print }' > \
      ${ORAL_BASE}_unmapped.alt
fi

# Scoring parameters found counting Alignment Score from bwakit hg38DH.fa.alt; this generates more supplementary alignments and missed odd MapQ 30 line
if [ "${VERSION_BASE}" = "hg38" ] && [ ! -e "additional_hg38_p14_${VERSION_HLA}_contigs.alt" ]; then
  if [ ! -e GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna.sa ]; then
    gzip -cd GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    bedtools maskfasta -fullHeader -fi GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -fo GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna -bed GCA_000001405.15_GRCh38_GRC_exclusions.bed
    bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna
  fi
  cat hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz GRCh38Patch14.fa.gz hla_gen.${VERSION_HLA}.fasta.gz > additional_hg38_p14_${VERSION_HLA}_contigs.fa.gz
  bwa mem -t`nproc` -A2 -B3 -O4 -E1 GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna additional_hg38_p14_${VERSION_HLA}_contigs.fa.gz \
    | samtools view -q60 - \
    | gawk '{ OFS="\t"; $10 = "*"; print }' > additional_hg38_p14_${VERSION_HLA}_contigs.alt
fi

# For the chm3 reference we don't yet have any defined alt sequences besides HLA; although we might consider new T2T references
if [ "${VERSION_BASE}" = "chm13" ] && [ ! -e "additional_chm13v2.0_${VERSION_HLA}_contigs.alt" ]; then
  if [ ! -e ${VERSION_BASE}${VERSION_PATCH}.fa.gz.sa ]; then
    bwa index chm13v2.0_maskedY_rCRS.fa.gz
  fi
  bwa mem -t`nproc` -A2 -B3 -O4 -E1 chm13v2.0_maskedY_rCRS.fa.gz hla_gen.${VERSION_HLA}.fasta.gz \
    | samtools view -q60 - \
    | gawk '{ OFS="\t"; $10 = "*"; print }' > additional_chm13v2.0_${VERSION_HLA}_contigs.alt
fi

if [ "${VERSION_BASE}" = "hg38" ] && [ ! -e "${VERSION}.fa.sa" ]; then
  cat GCA_000001405.15_GRCh38_full_analysis_set.fna.alt \
      ${DECOY_BASE}_unmapped.alt \
      additional_hg38_p14_${VERSION_HLA}_contigs.alt \
      ${ORAL_BASE}_unmapped.alt \
    > ${VERSION}.fa.alt

  zcat GCA_000001405.15_GRCh38_full_analysis_set_masked.fna.gz ${DECOY_BASE}_unmapped.fna.gz \
       hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz GRCh38Patch14.fa.gz hla_gen.$VERSION_HLA.fasta.gz ${ORAL_BASE}_unmapped.fna.gz > ${VERSION}.fa
  bwa index ${VERSION}.fa

  samtools faidx ${VERSION}.fa
  samtools dict -a "GRCh38" -s "Homo Sapiens" -u "${VERSION}.fa" ${VERSION}.fa -o ${VERSION}.dict
fi

if [ "${VERSION_BASE}" = "chm13" ] && [ ! -e "${VERSION}.fa.sa" ]; then
  cat ${DECOY_BASE}_unmapped.alt \
      additional_chm13v2.0_${VERSION_HLA}_contigs.alt \
      ${ORAL_BASE}_unmapped.alt \
    > ${VERSION}.fa.alt

  zcat chm13v2.0_maskedY_rCRS.fa.gz ${DECOY_BASE}_unmapped.fna.gz \
       hla_gen.$VERSION_HLA.fasta.gz ${ORAL_BASE}_unmapped.fna.gz > ${VERSION}.fa
  bwa index ${VERSION}.fa

  samtools faidx ${VERSION}.fa
  samtools dict -a "chm13v2.0" -s "Homo Sapiens" -u "${VERSION}.fa" ${VERSION}.fa -o ${VERSION}.dict
fi

#gatk-4.1.2.0/gatk FindBadGenomicKmersSpark -R ${VERSION}.fa -O ${VERSION}.fa.txt
