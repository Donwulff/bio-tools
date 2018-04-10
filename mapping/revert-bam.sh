#!/bin/sh

# Tested to work on Windows 10 WSL Ubuntu; install a java JRE
# Leaves nh and qf tags in BigY BAM; they seem properietary, so they shouldn't hurt

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
cores=8
coremem=$((javamem/$cores))

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

exit

# WARNING: This is NOT according to the Broad Institute GATK 4.0 Best Practices, but leaving it here for reference
(
samtools fastq -tn $outfile \
  | bwa mem -p -t $cores -M hs38DH.fa - \
  | samtools view -b -o ${bamfile%%.bam}.mem.bam
# Unfortunately, MarkDuplicates seeks back to beginning of the BAM so alignment can't just be piped in; improve later
java -jar picard.jar MarkDuplicates INPUT=${outfile%%.bam}.mem.bam OUTPUT=/dev/stdout METRICS_FILE=${bamfile}.dup \
    TAGGING_POLICY=All OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname COMPRESSION_LEVEL=0 \
    READ_NAME_REGEX="CL10.......L.C([0-9]+)R([0-9]+)_([0-9]+)" \
  | samtools sort -@$cores -m${coremem}G -l9 -o ${bamfile%%.bam}.srt.bam
) 2>&1 | tee ${bamfile}.log
