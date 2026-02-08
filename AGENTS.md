# Repository Guidelines

## Project Structure & Module Organization
- `annotate/`: VCF annotation and manipulation scripts (e.g., `annotate_y.sh`, `genos_annotate.sh`).
- `mapping/`: Core BAM processing workflows and reference build scripts (`revert-bam.sh`, `BQSR.sh`, `GRCh38_bwa_index.sh`) plus large reference assets and `index/` outputs.
- `util/`: Helper utilities and environment setup (`build_env.sh`, Perl helpers).
- `results/`: Example outputs (PDFs) for reference/comparison.

Scripts are standalone and typically run from their own directory. Large reference files live under `mapping/`; avoid editing or moving them unless required.

## Build, Test, and Development Commands
- `./util/build_env.sh`: Build common bioinformatics dependencies from source (samtools, bcftools, bwa, htslib, fastp, etc.).
- `cd mapping && ./GRCh38_bwa_index.sh`: Build reference genome indices (writes large outputs to `mapping/index/`).
- `./mapping/revert-bam.sh input.bam [optional_gigabytes_of_memory]`: Main pipeline to create uBAM, align, mark duplicates, and sort.
- `./mapping/BQSR.sh sample.sorted.bam [optional_gigabytes_of_memory]`: Run BQSR on a sorted BAM.

## Coding Style & Naming Conventions
- Shell scripts use POSIX `sh` (`#!/bin/sh`) and uppercase configuration variables (e.g., `REF`, `COMPRESS`).
- Perl utilities use `.pl` suffixes.
- No repo-wide formatter is enforced; match existing script style and prefer clear, inline comments only where needed.

## Testing Guidelines
- There is no centralized test runner in the repository.
- Validate changes by running the affected script on representative data and checking outputs with domain tools (e.g., `samtools`, `bcftools`).
- Use `results/` PDFs as qualitative comparisons where relevant.

## Commit & Pull Request Guidelines
- Existing history favors short, imperative summaries (sometimes with `fix:` prefixes). Keep the first line concise and scoped.
- PRs should include:
  - A brief description of the workflow impact.
  - Inputs/outputs or sample command lines used to validate.
  - Notes about large downloads, new reference assets, or configuration changes.

## Configuration & Environment Tips
- Scripts look for `./bio-tools.cfg` in the working directory for local overrides.
- Use `TMP=/path/to/ssd` to point temporary files to fast storage.
- Java 8+ is required for Picard/GATK steps.
