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
- There is no repository-wide test runner. `sh annotate/run_tests.sh` runs the focused annotation
  suites (synthetic comparator/reference/resource/backend/batch/resolver tests and the local-tree checks).
- Validate changes by running the affected script on representative data and checking outputs with domain tools (e.g., `samtools`, `bcftools`).
- Use `results/` PDFs as qualitative comparisons where relevant.
- Synthetic marker regressions: `python3 annotate/test_y_haplo_from_markers.py`; local-tree checks:
  `python3 annotate/test_ytree.py`. Neither requires opening private sample files.

## Commit & Pull Request Guidelines
- Standing user instruction: commit completed, validated work automatically in focused patches or
  batches. Leave unfinished work uncommitted and never mix unrelated user edits into a batch.
  Do not push or merge remote PRs unless explicitly requested.
- Existing history favors short, imperative summaries (sometimes with `fix:` prefixes). Keep the first line concise and scoped.
- PRs should include:
  - A brief description of the workflow impact.
  - Inputs/outputs or sample command lines used to validate.
  - Notes about large downloads, new reference assets, or configuration changes.

## Configuration & Environment Tips
- Liftover drivers, reference staging and marker/chain fetch helpers share a loader for trusted
  `./bio-tools.cfg` (or `BIO_TOOLS_CFG`). It is executable POSIX shell: load it before defaults/CLI
  parsing, so explicit options win. A missing explicit config fails; a missing default is allowed.
  Relative paths are analysis-working-directory relative, not config-file relative. Other workflows
  still have their own config handling; do not assume a repository-wide resource registry.
  Active `bio-tools.cfg` files are gitignored; commit only clearly named generic examples/templates.
- Key variables:
  - `Y_RESOURCES_DIR`: Preferred location for YBrowse marker files (snps_hg38.vcf.gz / .gff3) and related data; used by fetch and experiment scripts.
  - `BIO_TOOLS_INDEX_DIR` (or `INDEX_DIR`): Root for discovering reference FASTAs (with .fai/.dict) for liftover and mapping scripts. Falls back to `./resources/ref` or `mapping/index`.
  - `TMP=/path/to/ssd`: Point temporary files (Picard, sorting, etc.) to fast storage.
  - `LIFTOVER_BACKEND`: Shared liftover/dictionary backend: `auto`, `gatk`, `picard-jar`, `picard-cli`.
    Explicit CLI mode wins. Auto prefers a supplied `PICARD_JAR`, GATK, a discovered jar, then Picard CLI.
- Use the Java version required by the selected Picard/GATK release. Liftover drivers and reference
  preparation accept `--java-opts` or `PICARD_JAVA_OPTS` for heap tuning; an explicit empty value
  supplies no extra JVM options.

## Documentation
In addition to this file:
- `SITUATION.md` — dated situation review, branch/PR inventory, evidence precedence and open work;
  start here before reviving an old task from a chat or an early working note.
- Top-level `README.md` (high-level).
- `annotate/README.md`, `mapping/README.md`, `util/README.md` — detailed per-module usage and notes.
- `annotate/LIFTOVER.md` — reference preparation, Picard gVCF limitations and comparison caveats.
- `METHODS.md` — practical remapping methods (especially aDNA).
- `NOTES.md` — findings, gotchas, troubleshooting (pairing, MarkDuplicatesSpark, EAGER interaction, etc.).
- `RUNLOG.md` — run history.
- `TODO.md` — specific improvement items.
- `iceman-y/README.md` — the Ötzi Y verdict, and the only page written to be read on its own. `NOTES.md`
  and `RUNLOG.md` deliberately still carry the Iceman working notes and commands: they predate this
  investigation, they are shared with the mapping work, and there is no clean line to cut them on.
- `iceman-y/research/` — dated literature/candidate and research-priority assessments imported from
  reviewed Work PR revisions; planning inputs, not new sample results or verified scheduler state.
- `mapping/CLAUDE.md` — older Claude-oriented overview (complements this file).

## Y/Liftover Guardrails (Important)
- Separate competitive whole-genome alignment from a smaller liftover target/reference query.
  A chrY-only query is not authorization to remap reads against chrY alone.
- Prefer non-extended references for Y workflows unless explicitly testing decoy behavior:
  - hs1: `chm13v2.0_maskedY_rCRS.fa(.gz)`
  - hg38: `GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna`
- Extended references (`*H3630*`, `*DH3630*`, `*O1102*`, `hs38d1`, oral/metagenome additions) are valid for some mapping tasks but are easy to misuse for marker liftover experiments.
- Batch and dual-Y reference discovery share `annotate/liftover_resources.py`. Both require
  `--allow-extended-ref` for known extended filenames; this is not a full content/assembly check.
  Root order: configured index root, `$Y_RESOURCES_DIR/ref`, analysis `resources/ref`, repository
  `mapping/index`. Do not reintroduce whitespace-split roots or independent candidate lists.
  Staging defaults to the neutral `hg38.fa` name, which does not imply primary-only sequence content.
- Picard/GATK liftover requires both `.fai` and `.dict` sidecars.
- Backend selection/JVM/tool arguments share `annotate/liftover_backend.py`; do not reintroduce
  independent shell invocations or unquoted JVM-option expansion. GATK mode rejects whitespace
  inside JVM values; use jar mode for those. Picard CLI uses child-scoped `JAVA_TOOL_OPTIONS`, which
  wrappers may override. See `annotate/LIFTOVER.md` before claiming exact heap/version provenance.
- Sidecar checks share `annotate/reference_utils.py` across marker queries and shell liftover
  helpers. Invalid existing sidecars fail without replacement; `.gz` does not prove BGZF.
  Defaults check metadata and small probes, not whole-sequence identity. Use
  `annotate/check_reference.py --check-md5` for an explicit full sequence/layout check.
- Use a fresh staging directory for changed/copied/subset references. Unchanged source symlinks
  may be reused; copied byte-offset indexes must not be reused after a layout change.
- For compressed FASTA (`.fa.gz` / `.fna.gz`), keep dictionary naming compatible:
  - `<ref_without_fasta_suffix>.dict` (e.g. `x.fa.gz -> x.dict`)
  - and alias `<ref_without_last_suffix>.dict` when needed (e.g. `x.fa.dict`)
- For modern/private samples: keep run outputs outside repo and do not commit sample IDs or per-sample results.
- Do not treat reject-heavy Picard gVCF output as validated block conversion; see `annotate/LIFTOVER.md`.
  The ordinary-VCF batch driver rejects confidence records by content. Use a new/empty output
  directory; failed stages/logs are retained, and only validated complete batches are published.
- Comparator IDs across assemblies need a verified crosswalk. Complement permission is not
  identity verification; callable agreement is not accuracy and alias counts are not molecules.
- VCF label rankings are heuristic. For the Iceman verdict use committed read-level evidence and
  `markers/tree_local.tsv`; a derived SNP is not by itself a separate downstream node.
- Preserve dated corrections and preregistered outcomes. Record PR-only work as pending integration;
  do not silently turn metadata research, no-power outcomes or planned tests into completed results.

## Data Handling & Privacy
- Ancient/public datasets (for example, published Iceman data) may be documented in-repo and committed.
- Modern/private human datasets must stay out of git history.
- Never commit modern sample identifiers, read names, VCF/gVCF contents, or derived per-sample reports.
- In docs/scripts, use placeholders like `<sample_id>` and `<private_analysis_dir>` for modern workflows.
- Keep private run outputs in external analysis directories (outside this repository); commit only reusable tooling and generic method notes.
