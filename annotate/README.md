## annotate - VCF file annotation and manipulation

### annotate_y.sh
Historical script used for annotating Tyrolean Iceman GRCH38 mapped Y-chromosome with ISOGG gbrowse data.
Illustrates liftover from Hg19 to Hg38 for features from gbrowse exported annotations tracks.
Finally, creates tabix file from genomic feature file and uses it to annotate sample variants.
Has not been tested with current output files.

### genos_annotate.sh
Genos Research provided exome VCF's have weird format which cointains NT=Not Targeted Regions and NC=No Call.
This script is intended to merge the Genos Research VCF with the variants in dbSNP to fill in ref-calls.
Different bcftools have had different behavior for merging SNP's, so this needs to be tested.
Two different merges are attempted, one for multiallelic and one for separate calls.
After merging, output is annotated with CADD and MCAP pathogenity values.
This script has some additional bits to try to check concordance against DNA.Land imputed genomes to see if we got merging right.

### install_htslib.sh
Installs prerequisites, gets and runs Thomas Krahn's BigY2 (hg38) annotation script.
Tested on Windows 10 with Windows Store installed WSL Ubuntu.
This script assumes you have only one BigY VCF zip in your Windows download folder(s).

### combine_tests.sh

### yoruban_to_rcrs.sh - Yoruban YRI mtDNA reference VCF to rCRS calls converter
