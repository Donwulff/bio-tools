# bio-tools
Bioinformatics tools and snippets

Updated mapping/alignment reference genome creation script to GRCh38.p14. Release p14 posed some challenges, as it includes two new revisions of previously released contigs. This has been solved by including both versions of the contigs (you need to be using alt-aware aligner anyway), and using samtools faidx to pick and choose desired contigs. As previously, contigs are in the order of the first patch that introduced them; this way new patches are added at the end of existing ones, for those workflows which depend on consistent ordering.

As a reminder, the updated reference genome is intended for experimenting, and should not be used for real without through validation for intended purpose. It's based off the GRCh38DH reference genome (bwa.kit and 1000 Genomes Project Phase 3) with hs38d1 decoy sequences, updated with latest HLA sequences and reference genome patches up to p14 as alt-contigs. Optionally, it includes the Human Oral Microbiome Database sequences to use as decoys/analysing saliva sequences. Decoy and HOMD sequences are filtered, leaving out sequences with longer than 100bp exact matches to the patched human genome and HLA.

TODO list:
* Oral microbiome decoy and/or classification support for the saliva sequences.
* Make revert-bam.sh recognize unmapped bam from content and/or extension and skip unmapping.
* See about workflow creating unmapped bam on the fly without needing to store it.
* Add some of my scripting for creating unmapped bam from fastq files automatically.
* Try to reproduce decoy sequence creation from BAC/fosmid clones to include new ones.
* Do same for de novo assembled Personal Genomes Project sequences.
* Validate and compare results with GIAB and others.
* Publish my variant calling scripts.

## [annotate](annotate) - VCF file annotation and manipulation

## [mapping](mapping) - Tools for producing analysis-ready BAM files

## [util](util) - Random utilities
