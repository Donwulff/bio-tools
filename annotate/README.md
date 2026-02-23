## annotate - VCF file annotation and manipulation

### annotate_y.sh
Historical script used for annotating Tyrolean Iceman GRCH38 mapped Y-chromosome with ISOGG gbrowse data.
Illustrates liftover from Hg19 to Hg38 for features from gbrowse exported annotations tracks.
Finally, creates tabix file from genomic feature file and uses it to annotate sample variants.
Has not been tested with current output files.

### y_haplo_from_vcf.sh
Current Y-chromosome marker workflow for VCF input.
Uses the YBrowse hg38 marker VCF (`snps_hg38.vcf.gz`) directly, merges markers with sample chrY calls,
filters derived markers via ancestral allele (`INFO/AA`), and reports deepest haplogroup-chain candidates.
Recommended for current Iceman and modern re-calls.
Supports caller-specific filtering via `--site-filter-mode` / `--site-filter-expr` (includes `deepvariant` mode).

### y_haplo_from_markers.py
Python marker resolver for both VCF and gVCF inputs.
Streams marker rows against sample rows, resolves each marker to `derived|ancestral|nocall|ambiguous`,
and writes:
- `<prefix>.marker_status.tsv` (all marker statuses)
- `<prefix>.derived.tsv` (derived-only rows, optional `--min-gq/--min-dp`)
- `<prefix>.summary.txt`

### y_clade_consistency.py
Scores candidate clades against derived marker tables (`*.derived.tsv`) and reports
support/conflict metrics per dataset. Intended for side-by-side comparison across
published/legacy/new analysis branches.

### y_path_rank.py
Ranks clade/path candidates from `*.marker_status.tsv` with tunable scoring:
- derived support
- ancestral conflict
- no-call handling
- optional down-weighting for possible deamination transitions (`C>T`, `G>A`).

### run_iceman_y_compare.sh
Single-command driver for reproducible multi-branch Iceman chrY comparison:
- runs marker-state extraction (`y_haplo_from_markers.py`)
- runs path ranking (`y_path_rank.py`)
- writes combined reports (`key_markers.tsv`, `subtree_status.tsv`, `top_paths.tsv`)
Defaults target:
- auto-detected local files (`./iceman.vcf`, `./iceman.gvcf`, `./chrY_called_Iceman_tst_hg38.vcf.gz`)
- marker fallback: `annotate/data/snps_hg38.vcf.gz`
- optional overrides via env vars: `YCOMPARE_MARKERS`, `YCOMPARE_ICEMAN_VCF`, `YCOMPARE_ICEMAN_GVCF`, `YCOMPARE_SOLID_VCF`

### run_modern_y_experiment.sh
Driver for modern chrY experiments with assembly-safety checks.
It compares `chrY` contig lengths between sample VCF and marker VCF and stops on mismatch
unless explicitly overridden.

Practical use:
- For hs1/CHM13 calls + hg38 markers, do **not** run direct assignment.
- First liftover marker sites to hs1 (or use native GRCh38 calls), then run marker-state extraction.
- Default output path is under the caller's working directory: `./experiments/...`

### fetch_ybrowse_markers.sh
Fetches the YBrowse hg38 marker file (`snps_hg38.vcf.gz`) into a chosen path
(default `./resources/snps_hg38.vcf.gz`), with optional fallback to an
existing local mirror if network fetch fails.

### run_y_dual_liftover_experiment.sh
Two-branch experiment driver for modern Y analysis:
- branch A: liftover markers `hg38 -> hs1`, then call against native hs1 sample
- branch B: liftover sample `hs1 -> hg38`, then call against hg38 markers

Default output path is under the caller's working directory: `./experiments/...`.
Run it from a private working directory for sensitive sample data.

### liftover_to_hg38_batch.sh
Batch liftover helper for hs1/CHM13 sample VCF/gVCF files to GRCh38.
Designed for generating GRCh38 callsets for downstream tools (Exomiser/PharmCAT).
Optional `--primary-only` emits a canonical-contig subset
(`chr1-22,chrX,chrY,chrM`) per lifted sample.

### compare_vcf_runs.sh
General run-vs-run VCF comparator for the same sample/reference.
Useful for DeepVariant A/B checks such as small-model on/off.
Outputs:
- all-sites and PASS-only overlap counts (`summary.tsv`)
- shared-site GT concordance (`shared_gt_diff`)
- per-run private/shared VCFs
- optional per-position comparison (`--regions-file`) for marker checks.

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
