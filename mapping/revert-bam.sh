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
# Currently gets and runs MarkDuplicatesSpark multi-threaded instead of Picard MarkDuplicates flow.
GATK_SPARK=4.1.2.0
# Don't delete bwa-mem mapped, raw, unsorted BAM file. Useful if you intend to test different processing.
KEEP_TEMPORARY=True

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

if [ ! -z "${GATK_SPARK}" ] && [ ! -e gatk-${GATK_SPARK}/gatk ];
then
  if ! echo "${SANITIZE}" | grep -qi '^Strue$';
  then
    echo "You can't use GATK_SPARK MarkDuplicatesSpark with single end reads, leave SANITIZE on."
    exit
  fi
  wget -nc https://github.com/broadinstitute/gatk/releases/download/${GATK_SPARK}/gatk-${GATK_SPARK}.zip
  unzip gatk-${GATK_SPARK}.zip
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

  echo "############ Reverting ${SAMPLE} into ${UBAMFILE} with Picard Tools"
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

  # Identify Complete Genomics/BGI/MGI style read names for library complexity estimation
  try_regex='CL10.......L.C([0-9]+)R([0-9]+)_([0-9]+)'
  if samtools view ${UBAMFILE} | head -1 | cut -f1 | grep -q -E "${try_regex}";
  then
    echo "############   Identified BGI/MGI/Complete Genomics style read names"
    if [ -z "${GATK_SPARK}" ];
    then
      regex="READ_NAME_REGEX=${try_regex}"
    else
      regex="--read-name-regex ${try_regex}"
    fi
  fi

  # Uses .hdr file also as a flag of whether mapping finished, in case we restart
  if [ ! -e ${BMAPFILE} ] && [ ! -e ${UBAMFILE}.hdr ];
  then
    check_space ${BMAPFILE}

    samtools view -H ${UBAMFILE} > ${UBAMFILE}.hdr
    # FTDNA BigY doesn't have PL tag, so add it for downstream processing.
    if ! grep -q "^@RG.*PL:" ${UBAMFILE}.hdr;
    then
      echo "############   Adding missing platform identifiers"
      sed -i "s/^@RG\t.*/&\tPL:ILLUMINA/" ${UBAMFILE}.hdr
    fi

    echo "############ Using BWA MEM to align ${UBAMFILE} against ${REF} into ${BMAPFILE}"
    # Casava 1.8 header, observed s/[12]:[YN]:[0-9]*:[^\/]*\/[12]\t//
    # According to Wikipedia this can also end in barcode or sample-ID, so removing any.
    samtools fastq -t ${UBAMFILE} \
      | eval sed "s/[12]:[YN]:[0-9]*:[^[:space:]]*[[:space:]]//" \
      ${FILTER} \
      | bwa mem ${BWAOPT} -p -t ${cores} -M -C -H ${UBAMFILE}.hdr ${REF} - \
      | samtools view -b -o ${BMAPFILE}
    rm ${UBAMFILE}.hdr
  fi

  echo "############ Marking duplicates and chromosome order sorting ${UBAMFILE} into ${SORTFILE}"
  check_space ${SORTFILE}

  # Unfortunately, MarkDuplicates seeks back to beginning of the input BAM so mapping can't just be piped in
  if [ -z "${GATK_SPARK}" ];
  then
  java -jar picard.jar MarkDuplicates INPUT=${BMAPFILE} OUTPUT=/dev/stdout METRICS_FILE=${SAMPLE}.dup \
    ASSUME_SORT_ORDER=queryname TAGGING_POLICY=All COMPRESSION_LEVEL=0 TMP_DIR=$tmp \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ${regex} \
      | java -jar picard.jar FifoBuffer BUFFER_SIZE=2147483645 DEBUG_FREQUENCY=61 \
      | samtools sort -T $tmp/${SAMPLE##*/} -@${cores} -m${percoremem}G -l${COMPRESS} \
      | tee ${SORTFILE} \
      | samtools index -@${cores} - ${SORTFILE}.bai
  else
    # This took 8:20 vs. 6:17 on 4 cores, before running out of memory in index generation. Spark temp file space about 2X BAM, lots of IO.
    # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_spark_transforms_markduplicates_MarkDuplicatesSpark.php
    # Results are indentical, but duplication metrics are less detailed, PG isn't added, and output compression level can't be separately set.
    # --bam-partition-size 33554432 is maximum & default
    rm -rf ${tmp}/MDS
    mkdir -p ${tmp}/MDS
    gatk-${GATK_SPARK}/gatk --java-options "-Xmx${javamem}G -Dsamjdk.compression_level=${COMPRESS}" MarkDuplicatesSpark -I ${BMAPFILE} -O ${SORTFILE} -M ${SAMPLE}.dup \
      --duplicate-tagging-policy All --tmp-dir ${tmp}/MDS --output-shard-tmp-dir ${tmp}/MDS/${SAMPLE##./}.parts --optical-duplicate-pixel-distance 2500 ${regex} 2>&1 \
        | grep -Ev "INFO (Executor|NewHadoopRDD|ShuffleBlockFetcherIterator|SparkHadoopMapRedUtil|FileOutputCommitter):" --line-buffered
    # At least for me, the index file generation fails even with enough memory, so let's just generate it for now.
    samtools index -@${cores} ${SORTFILE}
  fi
  # Contrary to MarkDuplicatesSpark documentation, NM and MD tags look fine, but if you need UQ tags you need to run SetNmMdAndUqTags
  # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_SetNmMdAndUqTags.php

  # set -e and set -o pipefail require bash, so check the final file is intact instead
  if [ -z ${KEEP_TEMPORARY} ] && samtools quickcheck ${SORTFILE};
  then
    echo "############   Deleting aligner temporary output file"
    rm ${BMAPFILE}
  fi
fi
