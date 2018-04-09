#!/bin/sh
bamfile=${1:-"sample1.bam"}
outfile=${bamfile%%.bam}.u.bam

if [ ! -e $bamfile ];
then
  echo "Runs Picard Tools on Java to create Broad Institute uBAM from input BAM"
  echo "$0 [input.bam] [gigabytes of memory to use, if not all available memory]"
  exit
fi

if [ ! -e picard.jar ];
then
  wget https://github.com/broadinstitute/picard/releases/download/2.18.2/picard.jar
fi

# Different fs from data files, prefer SSD
tmp=${TMP:-"/tmp"}

# Sizes in bytes
inputsize=`du -D -b $bamfile | gawk '{print $1}'`
tempsize=$((inputsize+(inputsize/2)))
outfree=`LC_ALL=C df -B1 $bamfile | grep -v "^File" | gawk '{print $4}'`
tempfree=`LC_ALL=C df -B1 $tmp | grep -v "^File" | gawk '{print $4}'`

if [ $tempsize -gt $tempfree ];
then
  echo "Approximately 1.5X input size is required for temporary storage $tmp"
  echo "Run with TMP=/path/to/tmp to use different path"
  exit
fi

if [ $inputsize -gt $outfree ];
then
  echo "Output file $outfile probably won't fit on remaining space"
  exit
fi

totalmem=`LC_ALL=C free | grep -e "^Mem:" | gawk '{print $7}'`
# Allow 2 gigabytes for runtime
javamem=${2:-$((totalmem/1024/1024-2))}
# https://sourceforge.net/p/picard/wiki/Main_Page/#q-a-picard-program-that-sorts-its-output-sambam-file-is-taking-a-very-long-time-andor-running-out-of-memory-what-can-i-do
bamrecords=$((javamem*250000))

# Ref: https://gatkforums.broadinstitute.org/gatk/discussion/6484
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
    TMP_DIR=$tmp
