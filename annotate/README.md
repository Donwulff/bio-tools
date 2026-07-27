## annotate - VCF file annotation and manipulation

### annotate_y.sh
Historical script used for annotating Tyrolean Iceman GRCH38 mapped Y-chromosome with ISOGG gbrowse data.
Illustrates liftover from Hg19 to Hg38 for features from gbrowse exported annotations tracks.
Finally, creates tabix file from genomic feature file and uses it to annotate sample variants.
Has not been tested with current output files.

### y_haplo_from_vcf.sh
Legacy chrY marker workflow for VCF input.
Current development uses `y_haplo_from_markers.py` + `y_path_rank.py` as the canonical engine.
Keep `y_haplo_from_vcf.sh` mainly for backward compatibility with older runbooks.

### y_haplo_from_markers.py
Python marker resolver for VCF/gVCF sample inputs and marker definitions from either VCF or GFF3.
Streams marker rows against sample rows, resolves each marker to `derived|ancestral|nocall|ambiguous`,
and writes:
- `<prefix>.marker_status.tsv` (all marker statuses)
- `<prefix>.derived.tsv` (derived-only rows, optional `--min-gq/--min-dp`)
- `<prefix>.summary.txt`
Notes:
- For marker GFF3 input, `yfull_node` is preferred for `HG` when available (fallback: `ycc_haplogroup`).
- Marker/source ordering does not need to be pre-sorted.

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
- supports `--auto-clade` to infer top-level clade before ranking subpaths.
- filters noisy mixed-branch labels containing `&`/`|`.
- adds synthetic SNP labels from marker IDs when labels are coarse
  (e.g. `HG=I1` + `ID=Y47125` can contribute to `I-Y47125` ranking).
- reports simple support-strength summaries from derived marker rows:
  `derived_dp_sum`, `derived_dp_mean`, `derived_dp_max`, `derived_gq_mean`.
- `--tree markers/tree_local.tsv` adds a topology-aware placement beside the label
  ranking, written to `--tree-out` (default `<out>.tree.tsv`). The two are independent:
  the ranking reads HG/ISOGG labels, the placement reads only marker names against the
  tree, so they can disagree. When they do, that is the result.

### ytree.py
Loader, validator and scorer for the local Y tree in `markers/tree_local.tsv` — the one
machine-readable statement of which node sits above which. Marker lists and coordinate
files never carried topology, so `G-Z6219` above `G-L166` and `G-Z6499` below it lived
only in prose until this existed.

- `load_tree()` validates on read and raises rather than guessing: one root, no cycles,
  every parent known, every status from the fixed vocabulary, every marker on exactly
  one node, every node with at least one defining marker.
- `TreeScorer` calls each node from its **defining** markers only. `equivalent` markers
  are counted and reported but can never deepen a placement. `refuted` nodes are scored
  and printed but never appear on a path.
- `status_from_call()` translates `ylib.site_call()`'s vocabulary. One ancestral read
  stays a nocall, as that function intends; it is counted separately and surfaces as a
  `weak_ancestral` caveat on any derived call it sits under.
- `to_newick()` exports the tree for external placement tools (pathPhynder/phynder).
  Branch lengths are counts of defining markers, not time.

Run directly to check the file or export it:
```bash
python3 annotate/ytree.py --newick tree.nwk --markers-out markers.tsv
```

### y_tree_place.py
Places samples on the local tree from read-level pileup tables (`y_markers_pileup.py`
output), rather than from VCF-derived label tables. This is the path used by the recent
Iceman/Oberbipp/Aesch work.

```bash
python3 annotate/y_tree_place.py --pileup iceman-y/results/testC/testC_*.tsv \
    --out iceman-y/results/placement/aesch_muttenz_placement.tsv
```
- columns are found by name; files without them are skipped with a note, so a whole
  results directory can be globbed
- `--sample NAME` for single-sample tables that carry no `sample` column
- rows are deduplicated on `(sample, marker)`. The block files overlap on purpose —
  `L166_defining` and `yfull_L166_defining` share ten positions — so without this every
  derived and ancestral total would be inflated by the overlap
- writes a per-node table and a one-row-per-sample summary carrying the terminal node,
  its status, what the sample is excluded from, and any caveats

### test_ytree.py
`python3 annotate/test_ytree.py` — no pytest, no fixtures on disk. Checks the real tree
file loads, that eight classes of malformed tree are refused, and that the scorer still
reproduces the placements this project made from reads (Iceman at `G-L166` while derived
at the refuted `G-Z6208`; Swiss Neolithic men at `G-Z6219`; Sardinians at `G-PF3239`).

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
Status: compatibility helper. Prefer `run_y_dual_liftover_experiment.sh`.

Practical use:
- For hs1/CHM13 calls + hg38 markers, do **not** run direct assignment.
- First liftover marker sites to hs1 (or use native GRCh38 calls), then run marker-state extraction.
- Default output path is under the caller's working directory: `./experiments/...`
- Treat modern sample IDs and outputs as private: keep runs outside the repo and do not commit per-sample artifacts.

### fetch_ybrowse_markers.sh
Fetches the YBrowse hg38 marker VCF (`snps_hg38.vcf.gz`) into a chosen path
(default `./resources/snps_hg38.vcf.gz`), with optional fallback to an
existing local mirror if network fetch fails (`--fallback` or env `YMARKERS_FALLBACK`).
The fetched file is sanitized to remove invalid ALT definitions that can break
strict VCF readers.
By default it also fetches marker GFF3 to the sibling path (`snps_hg38.gff3`);
use `--no-gff3` to skip.
If `./bio-tools.cfg` defines `Y_RESOURCES_DIR`, default output is
`$Y_RESOURCES_DIR/snps_hg38.vcf.gz`.

### sanitize_marker_vcf.sh
Cleans marker VCF allele fields:
- removes ALT alleles equal to REF
- removes duplicate ALT alleles within a record
- drops records left with no ALT alleles

### fetch_liftover_chains.sh
Fetches `hg38ToHs1` and `hs1ToHg38` chain files to `./resources/chains/`
using UCSC mirrors with MARBL chain-name fallback.

### lib_liftover.sh
Shared shell helpers used by liftover drivers (`prepare_y_refs.sh`,
`run_y_dual_liftover_experiment.sh`, `liftover_to_hg38_batch.sh`):
- Picard auto-detection
- FASTA sidecar preparation (`.fai`, `.dict`)
- compressed-FASTA dictionary alias handling (`ref.dict` and `ref.fa.dict`)

### prepare_y_refs.sh
Stages hs1 and hg38 reference FASTAs into `./resources/ref/` and ensures
required sidecars exist (`.fai`, `.dict`). By default it symlinks source
FASTA files; use `--copy` to copy. Use `--primary-only` to build staged
FASTA files containing canonical contigs only (`chr1-22,chrX,chrY,chrM` by default).

### run_y_dual_liftover_experiment.sh
Two-branch experiment driver for modern Y analysis:
- branch A: liftover markers `hg38 -> hs1`, then call against native hs1 sample
- branch B: liftover sample `hs1 -> hg38`, then call against hg38 markers
- branch C (if GFF3 markers are available): rank hg38-lifted sample with hg38 marker GFF3 (`yfull_node` labels) for terminal-style naming
- branch D: path ranking with auto top-level clade (default) or manual `--clade-prefix`
- branch E: condensed `top_paths.tsv` summary (top positive-score rows per dataset)

Default output path is under the caller's working directory: `./experiments/...`.
Run it from a private working directory for sensitive sample data.
Supports `--java-opts` (or env `PICARD_JAVA_OPTS`) for LiftoverVcf heap sizing.
Uses `RECOVER_SWAPPED_REF_ALT=true` to retain loci where reference/alternate
alleles are swapped between assemblies.
Accepts `--markers-gff3-hg38` to force a specific GFF3 marker file.
References can be passed directly from your index build (full assemblies are fine);
the script ensures `.fai` and `.dict` sidecars exist before liftover.
If refs are omitted, the script auto-detects from:
- `--ref-root` or `BIO_TOOLS_INDEX_DIR`/`INDEX_DIR`
- `./resources/ref`
- `mapping/index`
and prefers no-alt/smaller analysis references before `...DH3630O1102.fa`.
By default the script rejects references that look like extended H3630/decoy sets;
use `--allow-extended-ref` only when you explicitly want that behavior.

Hard-won operational notes:
- For routine Y analysis, prefer non-extended refs (`chm13v2.0_maskedY_rCRS.fa(.gz)` and
  `GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna`) to avoid decoy-driven confusion.
- Picard liftover is strict about reference sidecars; ensure both `.fai` and `.dict` exist.
- With compressed FASTA (`.fa.gz`), some tool paths expect `ref.dict` and others `ref.fa.dict`.
  This script now normalizes/aliases both dictionary names to avoid false "missing dictionary" failures.
Canonical usage:
- Male sample QC + assignment cross-check: use this script.
- It is the maintained modern Y entrypoint.

### liftover_to_hg38_batch.sh
Batch liftover helper for hs1/CHM13 sample VCF/gVCF files to GRCh38.
Designed for generating GRCh38 callsets for downstream tools (Exomiser/PharmCAT).
Optional `--primary-only` emits a canonical-contig subset
(`chr1-22,chrX,chrY,chrM`) per lifted sample.
Supports `--java-opts` (or env `PICARD_JAVA_OPTS`) for LiftoverVcf heap sizing.
Uses `RECOVER_SWAPPED_REF_ALT=true` to retain loci where reference/alternate
alleles are swapped between assemblies.
Canonical usage:
- For non-Y pipelines (including female exome), use this script for hs1 -> hg38 conversion.
- `--ref-hg38` can be omitted if your index root is discoverable via `--ref-root` or `BIO_TOOLS_INDEX_DIR`.

### compare_vcf_runs.sh
General run-vs-run VCF comparator for the same sample/reference.
Useful for DeepVariant A/B checks such as small-model on/off.
Outputs:
- all-sites and PASS-only overlap counts (`summary.tsv`)
- shared-site GT concordance (`shared_gt_diff`)
- per-run private/shared VCFs
- optional per-position comparison (`--regions-file`) for marker checks.

### compare_marker_status_sets.py
Comparator for `*.marker_status.tsv` outputs across multiple datasets.
Useful for side-by-side checks such as:
- hs1 native vs hs1->hg38 liftover vs direct hg38 callset
- filtered marker subsets (e.g. terminal/trunk SNP IDs)
Outputs:
- `<prefix>.matrix.tsv` row-wise status/GT/DP/GQ/HG across datasets
- `<prefix>.summary.tsv` per-dataset status counts + pairwise agreement/flip counts

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
