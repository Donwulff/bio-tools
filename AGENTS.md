# Repository Guidelines

## Project Structure & Module Organization
- `annotate/`: VCF annotation, Y-chromosome marker analysis, and cross-assembly liftover tools. Legacy scripts (`annotate_y.sh`, `genos_annotate.sh`, `y_haplo_from_vcf.sh`) coexist with the current canonical Python engine (`y_haplo_from_markers.py`, `y_path_rank.py`, `y_clade_consistency.py`, `compare_marker_status_sets.py`) and driver scripts for reproducible Iceman comparisons and modern dual-liftover hs1↔hg38 experiments. Also includes fetch/sanitize helpers for YBrowse markers and Picard liftover chain/reference preparation.
- `mapping/`: Core BAM processing workflows and reference build scripts (`revert-bam.sh`, `BQSR.sh`, `GRCh38_bwa_index.sh`) plus large reference assets, additional alt-contig lists, patch FASTAs, and `index/` outputs. Contains `CLAUDE.md` (Claude-specific session notes, kept for the same reason as the root AGENTS.md).
- `DeDup/`: Reference copy of the third-party DeDup (Gradle/Java) tool for PCR duplicate removal on merged reads in aDNA/EAGER data. The entire source tree is ignored (see `.gitignore`); this repo only references build notes and (optionally) a local jar. Does not belong in this repo long-term. See the Donwulff/DeDup fork for the maintained source + Gradle wrapper.
- `util/`: Helper utilities and environment setup (`build_env.sh`, Perl helpers such as `bindex.pl` and `cigar-hist.pl`, Python tools such as `eager_repair_bam.py`, `is_stats.py`, and `samtools_if.sh`). Running the build populates large `util/src/` trees containing full source + build artifacts for samtools, bcftools, bwa, htslib, bedtools, etc.
- `resources/`: Fetched public reference data such as YBrowse marker sets (`snps_hg38.vcf.gz` + `.gff3`), liftover chains, and staged reference FASTAs (populated by scripts like `fetch_ybrowse_markers.sh`, `prepare_y_refs.sh`, `fetch_liftover_chains.sh`). Can be redirected via config.
- `notebooks/`: Jupyter analysis templates (e.g. `otzi_reanalysis_template.ipynb`).
- `iceman-y/`: The Ötzi Y-haplogroup investigation — the public case this tooling is developed against, kept together so the tooling repo does not carry one analysis's paperwork in its root. `README.md` is the verdict; `prereg/` holds the pre-registrations; `results/` holds every table the verdict cites. Nothing here is a dependency of `annotate/`, `mapping/` or `util/` — the arrow runs the other way.
- `markers/`: Marker name lists, site tables, and `tree_local.tsv` — the machine-readable local Y tree (node, parent, status, defining/equivalent markers). Read by `annotate/ytree.py`; see `annotate/README.md`.
- `results/`: Example outputs (PDFs for QC) for reference/comparison. Iceman analysis output lives under `iceman-y/results/`.

Scripts are standalone and typically run from their own (sub)directory. Large/generated content lives under `mapping/index/`, `resources/`, `util/src/`, and `__pycache__/` (see the root `.gitignore` for the full maintained list). The DeDup/ third-party tree is fully ignored. Avoid committing or editing generated/vendored content unless explicitly required.

## Build, Test, and Development Commands
- `./util/build_env.sh`: Build common bioinformatics dependencies from source (samtools, bcftools, bwa, htslib, fastp, etc.). This creates large trees under `util/src/`.
- `cd mapping && ./GRCh38_bwa_index.sh`: Build reference genome indices (writes large outputs to `mapping/index/`). Supports hg38 patches and CHM13/hs1.
- `./mapping/revert-bam.sh input.bam [optional_gigabytes_of_memory]`: Main pipeline to create uBAM, align, mark duplicates, and sort.
- `./mapping/BQSR.sh sample.sorted.bam [optional_gigabytes_of_memory]`: Run BQSR on a sorted BAM.
- DeDup (aDNA dedup for EAGER data): Use pre-built jar, bioconda, or build from the maintained fork (https://github.com/Donwulff/DeDup — includes Gradle wrapper and fixes). See NOTES.md. The source tree is ignored here.

See also `TODO.md` for pending improvements to build scripts.

## Coding Style & Naming Conventions
- Shell scripts primarily use POSIX `sh` (`#!/bin/sh`); a minority use `bash`. Use uppercase configuration variables (e.g., `REF`, `COMPRESS`).
- Python scripts (`.py`) are first-class and implement the modern Y-chromosome marker resolution, haplogroup path ranking, and comparison logic.
- Perl utilities use `.pl` suffixes.
- DeDup/ is a Gradle + Java project.
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
- Scripts look for `./bio-tools.cfg` (or path in `BIO_TOOLS_CFG`) in the working directory and source it for local overrides. Several annotate/ drivers also support environment variables for paths.
- Key variables:
  - `Y_RESOURCES_DIR`: Preferred location for YBrowse marker files (snps_hg38.vcf.gz / .gff3) and related data; used by fetch and experiment scripts.
  - `BIO_TOOLS_INDEX_DIR` (or `INDEX_DIR`): Root for discovering reference FASTAs (with .fai/.dict) for liftover and mapping scripts. Falls back to `./resources/ref` or `mapping/index`.
  - `TMP=/path/to/ssd`: Point temporary files (Picard, sorting, etc.) to fast storage.
- Java 8+ is required for Picard/GATK steps. Some liftover scripts accept `--java-opts` or `PICARD_JAVA_OPTS` for heap tuning.

## Documentation
In addition to this file:
- Top-level `README.md` (high-level).
- `annotate/README.md`, `mapping/README.md`, `util/README.md` — detailed per-module usage and notes.
- `METHODS.md` — practical remapping methods (especially aDNA).
- `NOTES.md` — findings, gotchas, troubleshooting (pairing, MarkDuplicatesSpark, EAGER interaction, etc.).
- `RUNLOG.md` — run history.
- `TODO.md` — specific improvement items.
- `iceman-y/README.md` — the Ötzi Y verdict, and the only page written to be read on its own. `NOTES.md`
  and `RUNLOG.md` deliberately still carry the Iceman working notes and commands: they predate this
  investigation, they are shared with the mapping work, and there is no clean line to cut them on.
- `mapping/CLAUDE.md` — older Claude-oriented overview (complements this file).

## Y/Liftover Guardrails (Important)
- Prefer non-extended references for Y workflows unless explicitly testing decoy behavior:
  - hs1: `chm13v2.0_maskedY_rCRS.fa(.gz)`
  - hg38: `GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna`
- Extended references (`*H3630*`, `*DH3630*`, `*O1102*`, `hs38d1`, oral/metagenome additions) are valid for some mapping tasks but are easy to misuse for marker liftover experiments.
- Picard/GATK liftover requires both `.fai` and `.dict` sidecars.
- For compressed FASTA (`.fa.gz` / `.fna.gz`), keep dictionary naming compatible:
  - `<ref_without_fasta_suffix>.dict` (e.g. `x.fa.gz -> x.dict`)
  - and alias `<ref_without_last_suffix>.dict` when needed (e.g. `x.fa.dict`)
- For modern/private samples: keep run outputs outside repo and do not commit sample IDs or per-sample results.

## Data Handling & Privacy
- Ancient/public datasets (for example, published Iceman data) may be documented in-repo and committed.
- Modern/private human datasets must stay out of git history.
- Never commit modern sample identifiers, read names, VCF/gVCF contents, or derived per-sample reports.
- In docs/scripts, use placeholders like `<sample_id>` and `<private_analysis_dir>` for modern workflows.
- Keep private run outputs in external analysis directories (outside this repository); commit only reusable tooling and generic method notes.
