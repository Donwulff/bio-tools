#!/bin/sh

# Tested to work on Windows 10 WSL Ubuntu & Ubuntu Xenial
# Leaves nh and qf tags in BigY BAM; they seem properietary
# If if bwa index reference exists (f.e. bwakit), simple mapping is performed
# Mapping requires significantly more storage space, currently not checked

# DEPENDENCIES: gawk; modern java; samtools (htslib) and bwa + index from bwakit if mapping

# Some BigY bam files have Casava 1.8 data in the bam read column; not sure what to do with this as it's in the original bam.
# This script preserves it in unmapped bam and filters it out before processing, but the extra space and Casava header can be cleaned up with:
# (samtools view -H sample1.bam; samtools view sample1.bam | awk '{ OFS="\t"; split($1,q," "); $1=q[1]; print }') | samtools view -b -o sample1.cleaned.bam

# File locations
SAMPLE=${1:-sample1.bam}
BASENAME=${SAMPLE%%.bam}
UBAMFILE=${BASENAME}.unmapped.bam
BMAPFILE=${BASENAME}.bwamem.bam
SORTFILE=${BASENAME}.sorted.bam
REF=hs38DH.fa
COMPRESS=5

# Extra options for BWA, ie. long reads with "-x pacbio" etc.
BWAOPT=

# Ref: https://gatkforums.broadinstitute.org/gatk/discussion/6484
# From page: Additionally, we invoke the SANITIZE option to remove reads that cause problems for certain tools, e.g. MarkIlluminaAdapters.
# Downstream tools will have problems with paired reads with missing mates, duplicated records, and records with mismatches in length of bases and qualities.
# Any paired reads file subset for a genomic interval requires sanitizing to remove reads with lost mates that align outside of the interval.
SANITIZE=true

# Trim sequencing adapters exactly using paired end read overlap; disableother trimming
# Reference: https://www.ncbi.nlm.nih.gov/pubmed/30423086
if which fastp > /dev/null;
then
  FILTER="| fastp -QLGp --stdin --stdout --interleaved_in -j ${BASENAME}.fastp.json -h ${BASENAME}.fastp.html -R \"fastp report on ${BASENAME}\""
fi

# Cores and memory
totalmem=`LC_ALL=C free | grep -e "^Mem:" | gawk '{print $7}'`
# Allow 2 gigabytes for runtime
javamem=${2:-$((totalmem/1024/1024-2))}
# https://sourceforge.net/p/picard/wiki/Main_Page/#q-a-picard-program-that-sorts-its-output-sambam-file-is-taking-a-very-long-time-andor-running-out-of-memory-what-can-i-do
bamrecords=$((javamem*250000))
cores=`nproc`
# From Java 6 update 18 max. heap is 1/4th of physical memory, so we can split 3/4th between cores for sorting.
percoremem=$((javamem*3/4/cores))

# If bio-tools.cfg exists where the script is run from, over-ride settings from it.
if [ -e ./bio-tools.cfg ];
then
  . ./bio-tools.cfg
fi

# Different fs from data files, prefer SSD
tmp=${TMP:-"/tmp"}

## Combining BigY hs37 and GRCh38, I removed reference and SANITIZE=true and ran:
# TMP=/mnt/SSD/ ./revert-bam.sh sample1BigY.bam
# TMP=/mnt/SSD/ ./revert-bam.sh sample1BigY2.bam
## Merging leaves reads from both files, apparently even with -c
# samtools merge -@4 -c -f sample1BigYmerge.bam sample1BigY.unmapped.bam sample1BigY2.unmapped.bam
## RevertSam expects reads in specific order, Picard Tools queryname order works.
# java -Xmx30G -jar picard.jar SortSam I=sample1BigYmerge.bam O=sample1BigYmerge.qsort.bam SORT_ORDER=queryname
## Added reference but not SANITIZE that would remove all duplicates, bwa undocumented removed duplicate read-id's:
# TMP=/mnt/SSD/ ./revert-bam.sh sample1BigYmerge.qsort.bam

if [ ! -e ${SAMPLE} ];
then
  echo "Runs Picard Tools on Java to create Broad Institute uBAM from input BAM"
  echo "$0 [input.bam] [gigabytes of memory to use, if not all available memory]"
  exit
fi

# Get a tested version of Picard Tools
if [ ! -e picard.jar ];
then
  wget https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar
fi

check_space ()
{
  # Figure out sizes in bytes
  inputsize=`du -D -b ${SAMPLE} | gawk '{print $1}'`
  outfree=`LC_ALL=C df -B1 ${SAMPLE%/*} | grep -v "^File" | gawk '{print $4}'`

  if [ ${inputsize} -gt ${outfree} ];
  then
    echo "Output file $1 (${inputsize}) probably won't fit on remaining space (${outfree})"
    exit
  fi

  # Figure out sizes in bytes - rule of thumb temporary is usually 1.5x times input
  tempsize=$((inputsize+(inputsize/2)))
  tempfree=`LC_ALL=C df -B1 ${tmp} | grep -v "^File" | gawk '{print $4}'`

  if [ ${tempsize} -gt ${tempfree} ];
  then
    echo "Approximately 1.5X input size is required for temporary storage $tmp"
    echo "Run with TMP=/path/to/tmp to use different path"
    exit
  fi
}

if [ ! -e ${UBAMFILE} ];
then
  check_space ${UBAMFILE}

  echo "Reverting $bamfile into ${UBAMFILE} with Picard Tools"
  java -Xmx${javamem}G -jar picard.jar RevertSam \
    I=${SAMPLE} \
    O=${UBAMFILE} \
    SANITIZE=${SANITIZE} \
    MAX_DISCARD_FRACTION=0.005 \
    ATTRIBUTE_TO_CLEAR=XT \
    ATTRIBUTE_TO_CLEAR=XN \
    ATTRIBUTE_TO_CLEAR=AS \
    ATTRIBUTE_TO_CLEAR=OC \
    ATTRIBUTE_TO_CLEAR=OP \
    ATTRIBUTE_TO_CLEAR=XS \
    ATTRIBUTE_TO_CLEAR=XA \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=true \
    REMOVE_DUPLICATE_INFORMATION=true \
    REMOVE_ALIGNMENT_INFORMATION=true \
    MAX_RECORDS_IN_RAM=$bamrecords \
    COMPRESSION_LEVEL=${COMPRESS} \
    TMP_DIR=$tmp
fi

# WARNING: This is NOT according to the Broad Institute GATK 4.0 Best Practices, but leaving it here for reference
# Does not perform quality control, adapter trimming or base quality score recalibration
if [ -e ${REF} ];
then
  echo "Mapping, marking duplicates and chromosome order sorting ${UBAMFILE} into ${SORTFILE}"

  # Identify Complete Genomics style read names for library complexity estimation
  try_regex='CL10.......L.C([0-9]+)R([0-9]+)_([0-9]+)'
  if samtools view ${UBAMFILE} | head -1 | cut -f1 | grep -q -E "${try_regex}";
  then
    regex="READ_NAME_REGEX=${try_regex}"
  fi

  # Uses .hdr file also as a flag of whether mapping finished, in case we restart
  if [ ! -e ${BMAPFILE} ] && [ ! -e ${UBAMFILE}.hdr ];
  then
    check_space ${BMAPFILE}

    samtools view -H ${UBAMFILE} > ${UBAMFILE}.hdr
    # FTDNA BigY doesn't have PL tag, so add it for downstream processing.
    if ! grep -q "^@RG.*PL:" ${UBAMFILE}.hdr;
    then
      sed -i "s/^@RG\t.*/&\tPL:ILLUMINA/" ${UBAMFILE}.hdr
    fi

    # Casava 1.8 header, observed s/[12]:[YN]:[0-9]*:[^\/]*\/[12]\t//
    # According to Wikipedia this can also end in barcode or sample-ID, so removing any.
    samtools fastq -t ${UBAMFILE} \
      | eval sed "s/[12]:[YN]:[0-9]*:[^[:space:]]*[[:space:]]//" \
      ${FILTER} \
      | bwa mem ${BWAOPT} -p -t ${cores} -M -C -H ${UBAMFILE}.hdr ${REF} - \
      | samtools view -b -o ${BMAPFILE}
    rm ${UBAMFILE}.hdr
  fi

  check_space ${SORTFILE}

  # Unfortunately, MarkDuplicates seeks back to beginning of the input BAM so mapping can't just be piped in
  java -jar picard.jar MarkDuplicates INPUT=${BMAPFILE} OUTPUT=/dev/stdout METRICS_FILE=${SAMPLE}.dup \
    ASSUME_SORT_ORDER=queryname TAGGING_POLICY=All COMPRESSION_LEVEL=0 TMP_DIR=$tmp \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ${regex} \
      | samtools sort -T $tmp/${SAMPLE##*/} -@${cores} -m${percoremem}G -l${COMPRESS} -o ${SORTFILE}
  # set -e and set -o pipefail require bash, so check the final file is intact instead
  if samtools quickcheck ${SORTFILE};
  then
    rm ${BMAPFILE}
  fi
fi
