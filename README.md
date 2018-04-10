# bio-tools
Bioinformatics tools and snippets

## annotate - VCF file annotation and manipulation

### annotate_y.sh
Historical script used for annotating Tyrolean Iceman GRCH38 mapped Y-chromosome with ISOGG gbrowse data.
Illustrates liftover from Hg19 to Hg38 for features from gbrowse exported annotations tracks.
Finally, creates tabix file from genomic feature file and uses it to annotate sample variants.
Has not been tested with current output files.

### genos_annotate.sh
Genos Research provided exome VCF's have weird format which cointains NT=Not Targeted Regions and NC=No Call.
This script is intended to merge thhe Genos Research VCF with the variants in dbSNP to fill in ref-calls.
Different bcftools have had different behavior for merging SNP's, so this needs to be tested.
Two different merges are attempted, one for multiallelic and one for separate calls.
After merging, output is annotated with CADD and MCAP pathogenity values.
This script has some additional bits to try to check concordance against DNA.Land imputed genomes to see if we got merging right.

### install_htslib.sh
Installs prerequisites, gets and runs Thomas Krahn's BigY2 (hg38) annotation script.
Tested on Windows 10 with Windows Store installed WSL Ubuntu.
This script assumes you have only one BigY VCF zip in your Windows download folder(s).

## mapping - BAM file mapping and processing

### revert_bam.sh
Runs Picard Tools on Java to create Broad Institute uBAM (Unmapped BAM) from input BAM.
Tested to work on Windows 10 WSL Ubuntu & Ubuntu Xenial.
Leaves nh and qf tags in BigY BAM; they seem properietary.
If if bwa index reference exists (f.e. bwakit), simple mapping is performed.
Mapping requires significantly more storage space, currently not checked.

#### Pros:
* GATK 4.0 Best Practices examples only support uBAM input
* Stores read metadata in standard way, in particular read groups
* RevertSam sanitizes reads and makes sure there are no errors
* Picard MarkDuplicates new mode requires uBAM sort order

#### Cons:
* Many quality control and adapter-trimming tools require FASTQ
* Picard tools isn't multithreaded, so sorting is very slow
* MarkDuplicates in queryname order can't be split by chromosome order
* Does not seem to store Illumina Casava 1.8 flags, this should be easy fix

From Picard MarkDuplicates documentation: "However, when the input is query-sorted (actually query-grouped), then unmapped
mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads."
When tested, the implementatiton actually required queryname order, crashing on query-grouped input.
You could sort the FASTQ files into queryname order, but this would take as much effort as RevertSam.
If marking all the duplicates is sgnificant for your workflow, uBAM seems currently only option, even if it's somewhat forced.
