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
# VERSION_HLA controls HLA sequence inclusion:
#   ""       = No HLA sequences
#   "H"      = Latest IPD-IMGT/HLA version with HLA ID naming (HLA00001) - safe for all pipelines
#   "A"      = Latest IPD-IMGT/HLA version with allele naming (HLA-A*01:01:01:01) - has : and * chars
#   "H3580"  = Specific version 3580 with HLA ID naming
#   "A3580"  = Specific version 3580 with allele naming
# HLA IDs (H) are persistent accession numbers; allele names (A) can change when phylogeny updates.
# Most pipelines (GATK etc.) choke on colons and asterisks, so HLA ID naming (H) is recommended.
# Delete Allele_status.txt to force re-fetch of latest version when using "A" or "H".
# VERSION_DECOY '' means not to include decoy sequences,
# VERSION_EXTRA is free-form used for local modifications
# WARNING! Oral microbiome is experimental. Version from 10.01 are huge, and require up to 64GB memory; 9.15 will fit in 32GB memory.
VERSION_BASE="hg38"
VERSION_PATCH="p14"
VERSION_DECOY="D"
VERSION_HLA="H"
VERSION_ORAL="9.15"
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
# VERSION_HLA modes:
#   ""           = No HLA sequences
#   "A" or "H"   = Latest version (A=accession naming, H=HLA allele naming)
#   "A1234"      = Specific version with accession naming
#   "H1234"      = Specific version with HLA allele naming
case "$VERSION_HLA" in
  "")
    # No HLA - will be handled by conditional includes below
    VERSION_HLA_CHR=""
    VERSION_HLA_NUM=""
    ;;
  A|H)
    # Latest version - fetch version number from IPD-IMGT/HLA
    # Uses -nc to preserve existing Allele_status.txt; delete it to force update
    VERSION_HLA_CHR="$VERSION_HLA"
    wget -nc https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allele_status.txt
    VERSION_HLA_NUM="$(grep version Allele_status.txt | tr -cd '[0-9]')"
    VERSION_HLA="${VERSION_HLA_CHR}${VERSION_HLA_NUM}"
    ;;
  A*|H*)
    # Specific version - parse the naming style and version number
    VERSION_HLA_CHR=$(printf "%.1s" "$VERSION_HLA")
    VERSION_HLA_NUM="${VERSION_HLA#?}"
    ;;
  *)
    echo "Invalid VERSION_HLA '$VERSION_HLA': must be empty, 'A', 'H', or specific version like 'A3580'"
    exit 1
    ;;
esac

# Helper variables for filename construction
# These use just the version number (e.g., "3590") so intermediate files can be shared between
# A and H naming conventions - the collision detection doesn't depend on how sequences are named.
if [ -n "$VERSION_HLA_NUM" ]; then
  HLA_VER_DOT_PREFIX=".${VERSION_HLA_NUM}"    # e.g., ".3590" for filenames like hg38p14.3590.fa.gz
  HLA_VER_USCORE_PREFIX="_${VERSION_HLA_NUM}" # e.g., "_3590" for filenames like ..._hg38p14_3590_genomic
else
  HLA_VER_DOT_PREFIX=""
  HLA_VER_USCORE_PREFIX=""
fi
# These use full VERSION_HLA (e.g., "A3590") for files where the naming convention matters,
# such as alt-contig files which contain sequence names with the chosen nomenclature.
if [ -n "$VERSION_HLA" ]; then
  HLA_NAMING_USCORE_PREFIX="_${VERSION_HLA}" # e.g., "_A3590" for alt contig filenames
  HLA_FASTA="hla_gen.${VERSION_HLA}.fasta.gz"
else
  HLA_NAMING_USCORE_PREFIX=""
  HLA_FASTA=""
fi

if [ -n "$VERSION_HLA" ] && [ ! -e "hla_gen.${VERSION_HLA}.fasta.gz" ]; then
  wget -nc https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta -O hla_gen.${VERSION_HLA_NUM}.fasta
  # Convert the HLA FASTA sequence names and compress it, no longer using bwa-kit HLA allele notation by default because : and * mess up most tools!
  if [ "$VERSION_HLA_CHR" = "A" ]; then
    sed "s/^>HLA:HLA[0-9]* />HLA-/" hla_gen.${VERSION_HLA_NUM}.fasta | gzip -c > hla_gen.${VERSION_HLA}.fasta.gz
  else
    sed "s/^>HLA:/>/" hla_gen.${VERSION_HLA_NUM}.fasta | gzip -c > hla_gen.${VERSION_HLA}.fasta.gz
  fi
fi

VERSION_ORAL_CODE=O$(echo "$VERSION_ORAL" | tr -d '.')
VERSION=${VERSION_BASE}${VERSION_PATCH}${VERSION_DECOY}${VERSION_HLA}${VERSION_ORAL_CODE}${VERSION_EXTRA}

# Index used to filter decoys/oral sequences by checking for collisions with human genome + HLA.
# Uses HLA version number only (not naming style) so it can be shared between A/H builds.
COLLISION_INDEX=${VERSION_BASE}${VERSION_PATCH}${HLA_VER_DOT_PREFIX}.fa.gz

# Construct the collision detection index
if [ "${VERSION_BASE}" = "hg38" ] && [ ! -e "${COLLISION_INDEX}.sa" ]; then
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
      $HLA_FASTA \
    > ${COLLISION_INDEX}
  bwa index ${COLLISION_INDEX}
fi

if [ "${VERSION_BASE}" = "chm13" ] && [ ! -e "${COLLISION_INDEX}.sa" ]; then
  wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/${VERSION_BASE}${VERSION_PATCH}.fa.gz
  # Use concatenated gzip's for speed because this is temporary index.
  cat ${VERSION_BASE}${VERSION_PATCH}.fa.gz \
      $HLA_FASTA \
    > ${COLLISION_INDEX}
  bwa index ${COLLISION_INDEX}
fi

## Steps to clean up decoy sequences
# Uses HLA_VER_USCORE_PREFIX so results can be shared between A/H naming conventions
DECOY_BASE=GCA_000786075.2_hs38d1_${VERSION_BASE}${VERSION_PATCH}${HLA_VER_USCORE_PREFIX}_genomic
if [ "${VERSION_DECOY}" != "" ] && [ ! -e "${DECOY_BASE}_unmapped.alt" ]; then
  # Filter out decoys which map to the current assembly for 101bp or more
  wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz
  bwa mem -t`nproc` -k101 ${COLLISION_INDEX} GCA_000786075.2_hs38d1_genomic.fna.gz > ${DECOY_BASE}.sam

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
# Uses HLA_VER_DOT_PREFIX so results can be shared between A/H naming conventions
ORAL_BASE=oral_microbiome_${VERSION_BASE}${VERSION_PATCH}${HLA_VER_DOT_PREFIX}_${VERSION_ORAL_CODE}_genomic
if [ "${VERSION_ORAL}" != "" ] && [ ! -e "${ORAL_BASE}_unmapped.alt" ]; then
  # Filter out decoys which map to the current assembly for 101bp or more
  wget -nc https://www.homd.org/ftp/genomes/PROKKA/V${VERSION_ORAL}/fsa/ALL_genomes.fsa -O oral_microbiome_${VERSION_ORAL_CODE}.fsa
  bwa mem -t`nproc` -k101 ${COLLISION_INDEX} oral_microbiome_${VERSION_ORAL_CODE}.fsa > ${ORAL_BASE}.sam
  samtools view -f0x4 ${ORAL_BASE}.sam | cut -f1 > \
    ${ORAL_BASE}_unmapped.list

  # Use the unmapped contigs to select matching decoy sequences from the analysis set
  samtools faidx oral_microbiome_${VERSION_ORAL_CODE}.fsa
  samtools faidx oral_microbiome_${VERSION_ORAL_CODE}.fsa -r ${ORAL_BASE}_unmapped.list | \
    bgzip -c > ${ORAL_BASE}_unmapped.fna.gz

  # Generate alt lines from alignment of remaining, unmapped decoys against the primary assy
  samtools view -f0x4 ${ORAL_BASE}.sam | \
    gawk -v OFS="\t" '{ $6 = length($10)"M"; $10 = "*"; $11 = "*"; NF=11; print }' > \
      ${ORAL_BASE}_unmapped.alt
fi

# Scoring parameters found counting Alignment Score from bwakit hg38DH.fa.alt; this generates more supplementary alignments and missed odd MapQ 30 line
# Alt file contains sequence names with naming convention, so uses full HLA_NAMING_USCORE_PREFIX
PATCH_HLA_ALT_BASE=additional_hg38_p14${HLA_NAMING_USCORE_PREFIX}_contigs
if [ "${VERSION_BASE}" = "hg38" ] && [ ! -e "${PATCH_HLA_ALT_BASE}.alt" ]; then
  if [ ! -e GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna.sa ]; then
    gzip -cd GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    bedtools maskfasta -fullHeader -fi GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -fo GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna -bed GCA_000001405.15_GRCh38_GRC_exclusions.bed
    bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna
  fi
  cat hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz GRCh38Patch14.fa.gz $HLA_FASTA > ${PATCH_HLA_ALT_BASE}.fa.gz
  bwa mem -t`nproc` -A2 -B3 -O4 -E1 GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna ${PATCH_HLA_ALT_BASE}.fa.gz \
    | samtools view -q60 - \
    | gawk '{ OFS="\t"; $10 = "*"; print }' > ${PATCH_HLA_ALT_BASE}.alt
fi

# For the chm13 reference we don't yet have any defined alt sequences besides HLA; although we might consider new T2T references
# Skip entirely if no HLA since chm13 has no patch contigs
CHM13_HLA_ALT_BASE=additional_chm13v2.0${HLA_NAMING_USCORE_PREFIX}_contigs
if [ "${VERSION_BASE}" = "chm13" ] && [ -n "$VERSION_HLA" ] && [ ! -e "${CHM13_HLA_ALT_BASE}.alt" ]; then
  if [ ! -e ${VERSION_BASE}${VERSION_PATCH}.fa.gz.sa ]; then
    bwa index chm13v2.0_maskedY_rCRS.fa.gz
  fi
  bwa mem -t`nproc` -A2 -B3 -O4 -E1 chm13v2.0_maskedY_rCRS.fa.gz $HLA_FASTA \
    | samtools view -q60 - \
    | gawk '{ OFS="\t"; $10 = "*"; print }' > ${CHM13_HLA_ALT_BASE}.alt
fi

if [ "${VERSION_BASE}" = "hg38" ] && [ ! -e "${VERSION}.fa.sa" ]; then
  cat GCA_000001405.15_GRCh38_full_analysis_set.fna.alt \
      ${DECOY_BASE}_unmapped.alt \
      ${PATCH_HLA_ALT_BASE}.alt \
      ${ORAL_BASE}_unmapped.alt \
    > ${VERSION}.fa.alt

  zcat GCA_000001405.15_GRCh38_full_analysis_set_masked.fna.gz ${DECOY_BASE}_unmapped.fna.gz \
       hg38Patch11.fa.gz GRCh38Patch12.fa.gz GRCh38Patch13.fa.gz GRCh38Patch14.fa.gz $HLA_FASTA ${ORAL_BASE}_unmapped.fna.gz > ${VERSION}.fa
  bwa index ${VERSION}.fa

  samtools faidx ${VERSION}.fa
  samtools dict -a "GRCh38" -s "Homo Sapiens" -u "${VERSION}.fa" ${VERSION}.fa -o ${VERSION}.dict
fi

if [ "${VERSION_BASE}" = "chm13" ] && [ ! -e "${VERSION}.fa.sa" ]; then
  cat ${DECOY_BASE}_unmapped.alt \
      ${VERSION_HLA:+${CHM13_HLA_ALT_BASE}.alt} \
      ${ORAL_BASE}_unmapped.alt \
    > ${VERSION}.fa.alt

  zcat chm13v2.0_maskedY_rCRS.fa.gz ${DECOY_BASE}_unmapped.fna.gz \
       $HLA_FASTA ${ORAL_BASE}_unmapped.fna.gz > ${VERSION}.fa
  bwa index ${VERSION}.fa

  samtools faidx ${VERSION}.fa
  samtools dict -a "chm13v2.0" -s "Homo Sapiens" -u "${VERSION}.fa" ${VERSION}.fa -o ${VERSION}.dict
fi

#gatk-4.1.2.0/gatk FindBadGenomicKmersSpark -R ${VERSION}.fa -O ${VERSION}.fa.txt
