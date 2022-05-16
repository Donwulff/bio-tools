# Tools for producing analysis-ready BAM files

Updated mapping/alignment reference genome creation script to 
GRCh38.p14. Release p14 posed some challenges, as it includes two new 
revisions of previously released contigs. This has been solved by 
including both versions of the contigs (you need to be using alt-aware 
aligner anyway), and using samtools faidx to pick and choose desired 
contigs. As previously, contigs are in the order of the first patch that 
introduced them; this way new patches are added at the end of existing 
ones, for those workflows which depend on consistent ordering.

As a reminder, the updated reference genome is intended for 
experimenting, and should not be used for real without through 
validation for intended purpose.

The reference genome is based off Heng Li's GRCh38DH reference genome 
(bwa-kit and 1000 Genomes Project Phase 3) with hs38d1 decoy sequences, 
updated with latest HLA sequences and reference genome patches up to p14 
as alt-contigs, and Genome Resource Consortium declared false sequences 
hard-masked. Optionally, it includes the Human Oral Microbiome Database 
sequences to use as decoys/analysing saliva sequences. Decoy and HOMD 
sequences are filtered, leaving out sequences with longer than 100bp 
exact matches to the patched human genome and HLA.

TODO list:
 * Oral microbiome classification support for the saliva sequences.
 * Make revert-bam.sh recognize unmapped bam from content and/or extension and skip unmapping.
 * See about workflow creating unmapped bam on the fly without needing to store it.
 * Add levioSAM2 lift-over into workflows.
 * Test and integrate DRAGENMAP into the example scripts.
 * Try to reproduce decoy sequence creation from BAC/fosmid clones to include new ones.
 * Do same for de novo assembled Personal Genomes Project sequences.
 * Validate and compare results with GIAB and others.
 * Publish my variant calling example scripts.

## revert_bam.sh
Runs Picard Tools on Java to create Broad Institute uBAM (Unmapped BAM) from input BAM.
Tested to work on Windows 10 WSL Ubuntu & Ubuntu Xenial Xerus, Bionic Beaver, Focal Fossa and Jammy Jellyfish at times.
If bwa index reference exists (f.e. bwakit), (not so) simple mapping is performed, if available with:
- [fastp](https://github.com/OpenGene/fastp) is used for fast, exact (from overlap) adapter trimming.
- [bwa mem](https://github.com/lh3/bwa) 0.7.17 (r1188) is standard, but latest dev works fine too.
- [k8 JavaScript interpreter](https://github.com/attractivechaos/k8/releases) For running the bwa-postalt.js script, pre-compiled easiest.
- [bwa-postalt.js](https://github.com/lh3/bwa/tree/master/bwakit) On the bwa/bwa-kit repository, for post-processing alt-alignments.
- [Picard](https://broadinstitute.github.io/picard/) or [GATK](https://github.com/broadinstitute/gatk/releases) For marking duplicate reads
- [htslib](https://github.com/samtools/htslib) For compressing and indexing genomic data.
- [samtools](https://github.com/samtools/samtools) For processing & indexing Sequence Alignment Maps.

This began as an example of how to run RevertSam. It still mostly is an example of how to run RevertSam.

But that necessitated showing how to map uBAM with BWA MEM, and that in turn required handling duplicates and so forth...
If you require cluster support, GATK MarkDuplicateSpark is now supported by the script, but for full GATK Best Practices compliance look instead at
https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
That said, I believe these scripts to create GATK Best Practices compatible BAM's in very efficient manner.

NB. GATK Best Practices MergeBamAlignment which this script omits for speed has some unheralded filtering by default.
MAX_INSERTIONS_OR_DELETIONS in particular is set to 1 and will filter 
any reads with more than single insertion/deletion. This script does not do that, but I believe this difference to be beneficial.
UNMAP_CONTAMINANT_READS might be beneficial for saliva sample, but it's not set by default, and not done by this script.
MergeBamAlignment will copy the alignments back into the Unmapped BAM in the exact (read-name) order they originally were, regardless 
of the order the aligner outputs them. This may cause differences in some special cases, in particular when comparing two different
BAM's, since genomic coordinate sorting will preserve the existing order for reads starting at the same genetic coordinate.
https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment

### uBAM Pros:
* GATK 4.0 Best Practices examples only support uBAM input
* Stores read metadata in standard way, in particular read groups
* RevertSam sanitizes reads and makes sure there are no errors
* Picard MarkDuplicates new mode requires uBAM sort order

### uBAM Cons:
* Many quality control and adapter-trimming tools require FASTQ
* Picard tools isn't multithreaded, so sorting is very slow
* MarkDuplicates in queryname order can't be parallel-split by chromosome order
* Does not seem to store Illumina Casava 1.8 flags, this should be easy fix

### Original motivation - new queryname order MarkDuplicates workflow for full duplicate identification
From Picard MarkDuplicates documentation: "However, when the input is query-sorted (actually query-grouped), then unmapped
mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads."
When tested, the implementatiton actually required queryname order, crashing on query-grouped input.
You could sort the FASTQ files into queryname order, but this would take as much effort as full RevertSam.
If marking all the duplicates is sgnificant for your workflow, uBAM seems currently best option, even if it's somewhat forced
by Broad Institute.

### Automatically handles many situations and issues found in personal genome sequences
* Automatically detects and removes adapter read-through on any PE sequencing (Dante Labs)
* Missing platform tag in read group (BigY)
* CASAVA headers in read names (BigY)
* Missing read pairs in chromosome-filtered BAM (BigY)
* Complete Genomics style read-name handling (Dante Labs, Genos Research)

### Other features
* Attempts to check and warn about obviously insufficient disk space
* Automatically uses CPU threads & memory available
* Temporary files can be directed to a fast disk
* MarkDuplicates or MarkDuplicatesSpark depending on needs
* Configure output compression level according to BAM file intent
* Stages are piped without unneccessary I/O where possible
* Input BAM headers and tags are lifted over directly during alignment
* Uses samtools for fast parallel compression and indexing

### Performance comparisons

http://gigadb.org/dataset/100274

| System         | Zlib     | Alignment| MarkDuplicatesSpark | MarkDuplicates |
|----------------|----------|---------:|--------------------:|---------------:|
| AMD 4x 3.1 Ghz | Standard |  44:00 h |              8:20 h |         6:17 h |

## BQSR.sh - Base Quality Score Recalibration
Broad Institute BQSR: https://software.broadinstitute.org/gatk/documentation/article?id=11081

Sequencing machines generate accuracy likelihood for each nucleotide base sequenced, called Base Quality Score.
Each sequencing machine, and flowcell (location on the sequencing machine) has unique profile of errors depending on context.
BQSR excludes lists of common genetic variation and creates an context-dependent error profile for the sequence.
It then uses this empirically created error profile to calibrate all the base quality scores of each read.
Most sequencing comapanies return raw BAM from before the BQSR is ran, so you may not need to run it for third party analysis.

BQSR may not always be beneficial: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048557/

Analysis is done split to chromosomes, but writing output is yet done sequentially to avoid I/O of concatenating chromosomes.

### Known variant sites that are filtered out
https://software.broadinstitute.org/gatk/documentation/article.php?id=1247

Note you need different known variant lists for different reference genome builds. I'm working on hg38 only currently.

Due to interested in Y chromosome phylogeny (https://www.yfull.com/ etc.) this script grabs International Society of Genetic 
Genealogy https://isogg.org/ YBrowse https://ybrowse.org Y variant list by YSEQ https://www.yseq.net/ to augment known sites 
for Y chromosome. For Y-targeted sequencing like BigY https://www.familytreedna.com/products/y-dna or Y Elite https://www.fullgenomes.com/
Y chromosome is also our ONLY source of sequencing error profile.

In this BQSR report, black dots indicate the error covariates for the whole genome, or the error profile that is.
The red dots indicate the same error covariates for X chromosome (Haploid male chromosome). The blue dots, which indicate
error covariates for Y chromosome are something else entirely.

https://github.com/Donwulff/bio-tools/blob/master/results/DanteLabs_PE100.sorted.bam.YvsX.pdf

The Y chromosome is difficult to sequence due to being largely heterochromatin and repeats. From the technical principles of
BQSR there is reason to think that the error profile over the rest of the genome is the correct one (The genome sequencer,
after all, does not know it's sequencing Y chromosome and the technical underpinnings should be the same). This means that 
BQSR likely over-estimates the error over Y chromosome, and for Y-targeted sequencing like Big Y and Y Elite which don't 
have other chromosomes the impact may be significant. Since women do not possess Y chromosome, it may also lead to sex-biased
error profiles. You need to evaluate whether to include Y chromosome in the error profile. By default I've included it, because 
it may smooth out error due to some covariates in BigY/Y Elite.

Difficulty of sequencing Y chromosome: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC138936/

### Some performance statistics on WGS run:
AMD 4x 3.1Ghz: 4:45 h to construct model, 8:50 h to apply it single-threaded (Standard zlib)
Preserve original qualities: 164GiB, otherwise 103GiB

### Accidental covariates?
I had to interrupt the alignment run, and only after running BQSR remembered samtools merge requires -c -p parameters or it 
will create new read group & program groups for collisions. However, this shows the latter 2/3rds of the second read group 
has empirically much higher error rate than the first third. Flowcell runtime as a covariate? Maybe, but note that the
EstimatedQReported is somewhat down too so the sequencer already knows it's less reliable.

https://github.com/Donwulff/bio-tools/blob/master/results/BGISEQ-500_PE100.sorted.bam.pdf

Here red and blue are error covariates of Chromosome 20 only (To save time) before and after applying the whole-genome BQSR.
Note the blue dots are not quite flat, so significant residual error not explained by the existing covariates still exists.

Insert size is another self-evident covariant, but a recent paper put it down to science:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6393434/

Unfortunately GATK doesn't currently support these covariates, but I'd like to see how much of error they explain.

## GRCh38_bwa_index.sh - EXPERIMENTAL reference genome build script

Now with both GRCh38 patch 14 and Telomere-to-Telomere CHM13v2 support.
