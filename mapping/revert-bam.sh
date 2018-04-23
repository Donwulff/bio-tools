#!/bin/sh

# Tested to work on Windows 10 WSL Ubuntu & Ubuntu Xenial
# Leaves nh and qf tags in BigY BAM; they seem properietary
# If if bwa index reference exists (f.e. bwakit), simple mapping is performed
# Mapping requires significantly more storage space, currently not checked

# DEPENDENCIES: gawk; modern java; samtools (htslib) and bwa + index from bwakit if mapping

# Different fs from data files, prefer SSD
tmp=${TMP:-"/tmp"}

# File locations
bamfile=${1:-sample1.bam}
outfile=${bamfile%%.bam}.unmapped.bam
reference=hs38DH.fa
compress=5

if [ ! -e $bamfile ];
then
  echo "Runs Picard Tools on Java to create Broad Institute uBAM from input BAM"
  echo "$0 [input.bam] [gigabytes of memory to use, if not all available memory]"
  exit
fi

# Get a tested version of Picard Tools
if [ ! -e picard.jar ];
then
  wget https://github.com/broadinstitute/picard/releases/download/2.18.2/picard.jar
fi

check_space ()
{
  # Figure out sizes in bytes
  inputsize=`du -D -b $bamfile | gawk '{print $1}'`
  outfree=`LC_ALL=C df -B1 $bamfile | grep -v "^File" | gawk '{print $4}'`

  if [ $inputsize -gt $outfree ];
  then
    echo "Output file $1 probably won't fit on remaining space"
    exit
  fi

  # Figure out sizes in bytes - rule of thumb temporary is usually 1.5x times input
  tempsize=$((inputsize+(inputsize/2)))
  tempfree=`LC_ALL=C df -B1 $tmp | grep -v "^File" | gawk '{print $4}'`

  if [ $tempsize -gt $tempfree ];
  then
    echo "Approximately 1.5X input size is required for temporary storage $tmp"
    echo "Run with TMP=/path/to/tmp to use different path"
    exit
  fi
}

check_space $outfile

# Cores and memory
totalmem=`LC_ALL=C free | grep -e "^Mem:" | gawk '{print $7}'`
# Allow 2 gigabytes for runtime
javamem=${2:-$((totalmem/1024/1024-2))}
# https://sourceforge.net/p/picard/wiki/Main_Page/#q-a-picard-program-that-sorts-its-output-sambam-file-is-taking-a-very-long-time-andor-running-out-of-memory-what-can-i-do
bamrecords=$((javamem*250000))
cores=`nproc`
# From Java 6 update 18 max. heap is 1/4th of physical memory, so we can split 3/4th between cores for sorting.
percoremem=$((javamem*3/4/cores))

# Ref: https://gatkforums.broadinstitute.org/gatk/discussion/6484
if [ ! -e $outfile ];
then
  echo "Reverting $bamfile into $outfile with Picard Tools"
  java -Xmx${javamem}G -jar picard.jar RevertSam \
    I=$bamfile \
    O=$outfile \
    SANITIZE=true \
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
    COMPRESSION_LEVEL=$compress \
    TMP_DIR=$tmp
fi

# WARNING: This is NOT according to the Broad Institute GATK 4.0 Best Practices, but leaving it here for reference
# Does not perform quality control, adapter trimming or base quality score recalibration
if [ -e $reference ];
then
  check_space ${bamfile%%.bam}.mem.bam

  echo "Mapping, marking duplicates and chromosome order sorting $outfile into ${bamfile%%.bam}.srt.bam"

  # Identify Complete Genomics style read names for library complexity estimation
  try_regex='CL10.......L.C([0-9]+)R([0-9]+)_([0-9]+)'
  if samtools view $outfile | head -1 | cut -f1 | grep -q -E "$try_regex";
  then
    regex="READ_NAME_REGEX=$try_regex"
  fi

  # Uses .hdr file also as a flag of whether mapping finished, in case we restart
  if [ ! -e ${bamfile%%.bam}.mem.bam ] && [ ! -e ${outfile}.hdr ];
  then
    samtools view -H $outfile > ${outfile}.hdr
    samtools fastq -t $outfile \
      | bwa mem -p -t $cores -M -C -H ${outfile}.hdr $reference - \
      | samtools view -b -o ${bamfile%%.bam}.mem.bam
    rm ${outfile}.hdr
  fi

  check_space ${bamfile%%.bam}.srt.bam

  # Unfortunately, MarkDuplicates seeks back to beginning of the input BAM so mapping can't just be piped in
  java -jar picard.jar MarkDuplicates INPUT=${bamfile%%.bam}.mem.bam OUTPUT=/dev/stdout METRICS_FILE=${bamfile}.dup \
    ASSUME_SORT_ORDER=queryname TAGGING_POLICY=All COMPRESSION_LEVEL=0 TMP_DIR=$tmp \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 $regex \
      | samtools sort -T $tmp/${bamfile##*/} -@$cores -m${percoremem}G -l${compress} -o ${bamfile%%.bam}.srt.bam
  # set -e and set -o pipefail require bash, so check the final file is intact instead
  if samtools quickcheck ${bamfile%%.bam}.srt.bam;
  then
    rm ${bamfile%%.bam}.mem.bam
  fi
fi
