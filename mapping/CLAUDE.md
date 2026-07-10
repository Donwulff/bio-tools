# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

bio-tools is a shell-based bioinformatics toolkit for creating GATK Best Practices-compliant, analysis-ready BAM files from next-generation sequencing data. The project supports multiple reference genomes (GRCh38, GRCh37, CHM13v2) and various sequencing platforms (Illumina, BGI/MGI, Complete Genomics).

## Repository Structure

- **mapping/** - Core BAM file creation and reference genome building (current directory)
- **annotate/** - VCF file annotation tools
- **util/** - Build environment setup and analysis utilities

## Key Scripts in mapping/

### revert-bam.sh
Main workflow: converts mapped BAM → unmapped BAM (uBAM) → BWA MEM alignment → duplicate marking → sorted output.
```bash
./revert-bam.sh input.bam [optional_gigabytes_of_memory]
```

### BQSR.sh
Base Quality Score Recalibration using GATK. Calibrates base quality scores using known variant sites.
```bash
./BQSR.sh sample.sorted.bam [optional_gigabytes_of_memory]
```

### GRCh38_bwa_index.sh
EXPERIMENTAL: Builds custom reference genomes with HLA sequences, alt-contigs, and optional HOMD decoys. Configuration via variables at top of script:
- `VERSION_BASE` - "hg38" or "chm13"
- `VERSION_PATCH` - patch version (e.g., "p14")
- `VERSION_HLA` - "A" for accession nomenclature
- `VERSION_ORAL` - HOMD version or "" to disable

## Build/Development Commands

**Setup development environment** (builds samtools, bcftools, bwa, htslib, fastp, etc. from source):
```bash
./util/build_env.sh
```

**Build reference genome index** (~179GB output in index/):
```bash
cd mapping && ./GRCh38_bwa_index.sh
```

## Configuration

Scripts look for `./bio-tools.cfg` in the working directory for local configuration overrides. Scripts auto-detect available CPU cores and memory. Environment variable `TMP=/path/to/ssd` can direct temporary files to fast storage.

## Architecture Notes

- Each script is standalone and independently executable
- Progressive pipeline: revert → align → deduplicate → sort → (optional) BQSR
- Supports both Picard and GATK Spark backends (Spark for multi-threading)
- Scripts automatically download required tools (Picard, GATK, etc.) on first run
- Java 1.8+ required for Picard/GATK operations

## Important Technical Details

- The scripts deliberately skip MergeBamAlignment (which filters reads with MAX_INSERTIONS_OR_DELETIONS=1)
- Alt-contig aware alignment requires bwa-postalt.js with k8 interpreter
- Y chromosome handling uses YBrowse variant list to augment known sites for BQSR
- Reference genome maintains consistent contig ordering across patches (new contigs appended at end)
