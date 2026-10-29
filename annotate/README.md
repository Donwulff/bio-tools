## annotate - VCF file annotation and manipulation

See [LIFTOVER.md](LIFTOVER.md) for the reference/conversion guide and known comparison limits, and
[the situation review](../SITUATION.md) for evidence status and unmerged Work PRs. Current Iceman
placement uses read-level pileups and the local tree. VCF label rankings remain exploratory.

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
Resolves marker rows to `derived|ancestral|nocall|ambiguous|filtered|other|unknown`,
and writes:
- `<prefix>.marker_status.tsv` (all marker statuses)
- `<prefix>.derived.tsv` (derived-only rows after filtering)
- `<prefix>.summary.txt`
Notes:
- For marker GFF3 input, `yfull_node` is preferred for `HG` when available (fallback: `ycc_haplogroup`).
- Marker/source ordering does not need to be pre-sorted.
- `--reference` supplies the sample assembly's indexed FASTA for interior gVCF block bases; it
  requires pysam and existing nonempty `.fai` (plus `.gzi` for BGZF). No reference index is created
  implicitly. Both `<*>` and `<NON_REF>` reference blocks are supported. Without a reference,
  interior block markers are `nocall` with reason `reference_required`, not copies of the anchor REF.
- `--min-gq/--min-dp` now apply before **all** output states and summaries. Missing or sub-threshold
  values and excluded site filters produce `filtered`, retaining GT/GQ/DP and the exclusion reason.
  The ranker counts these as no-calls (default weight zero); the comparator retains `filtered`.
- The original 13 columns remain in order; column 14 is `reason` (`.` for a resolved call).
  `marker_status.tsv` and `derived.tsv` remain headerless. Older 13-column tables remain readable by
  the ranker/comparator, but contain no new filter/reference validation and must not be relabelled.
- For blocks, DP reports conservative `MIN_DP` when present, otherwise DP. Zero effective block
  depth is not allele evidence; indel markers in blocks remain explicit no-calls. DeepVariant mode
  permits `FILTER=.` on reference blocks, but not on ordinary variant records.
- No sample analyses were rerun for this correction. Use a fresh output prefix for new results;
  see [LIFTOVER.md](LIFTOVER.md#comparing-native-hs1-lifted-hg38-and-direct-hg38).

Focused annotation regression checks (no sample inputs; indexed-reference checks need pysam):

```bash
sh annotate/run_tests.sh
```

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
- These DP statistics are depth at derived marker rows, not allele-specific read support or counts
  of independent molecules; the resolver does not export AD. Scores and `support_ratio` are
  descriptive heuristics, not probabilities. Equivalent names can count the same locus repeatedly.
- Candidate labels synthesized from IDs help locate individual SNPs; they do not establish that
  each SNP defines a distinct branch. Deamination weighting defaults to 0.5 even for modern data.
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
or `$Y_RESOURCES_DIR/chains/`, using UCSC mirrors with MARBL chain-name fallback.

### lib_liftover.sh
Shared shell helpers used by liftover drivers (`prepare_y_refs.sh`,
`run_y_dual_liftover_experiment.sh`, `liftover_to_hg38_batch.sh`):
- Picard auto-detection
- FASTA sidecar preparation (`.fai`, `.dict`)
- compressed-FASTA dictionary alias handling (`ref.dict` and `ref.fa.dict`)

Reference discovery is shared through `liftover_resources.py`; Picard backend discovery/execution
still needs consolidation. The shell helpers also load trusted local config before CLI parsing.
Sidecar validation now uses the same `reference_utils.py` implementation as marker queries. Empty,
malformed or mismatching existing sidecars fail rather than being overwritten. Set `SCRIPT_DIR`
to the annotate directory before sourcing the library. See [LIFTOVER.md](LIFTOVER.md).

### liftover_resources.py
One selector for the batch and dual-Y drivers: explicit reference paths, configured index/resource
roots, analysis-local `resources/ref`, then repository `mapping/index`. Paths are not whitespace-split
and no recursive filesystem search or download is performed by discovery. Base analysis references
are preferred; extended filenames require `--allow-extended-ref` in either driver. This is a filename
safeguard, not verification of the assembly or every contig. See
[resource configuration](LIFTOVER.md#resource-configuration) for precedence, defaults and limitations.

Both drivers, reference staging and the two fetch helpers now share the `BIO_TOOLS_CFG` loader.
An explicitly named missing config is an error; absence of default `./bio-tools.cfg` is allowed.
The file is executable POSIX shell, not an inert settings file: use only trusted local configs.

### check_reference.py / reference_utils.py
Read-only checks for FASTA format, FAI structure, BGZF/GZI compatibility, dictionary names/lengths/
order, dictionary-alias consistency and indexed boundary access. Requires pysam for indexed probes;
plain gzip is rejected, regardless of the filename suffix. The shell helper creates missing
sidecars with existing backends, then validates them; this checker never rebuilds an invalid one.

```bash
python3 annotate/check_reference.py /path/to/native.fa
python3 annotate/check_reference.py /path/to/native.fa --check-md5
```

Default checks do not prove whole-sequence identity. `--check-md5` explicitly scans every sequence
to verify dictionary M5 values and FAI byte layout, so it costs a full reference read. GZI validation
remains structural plus access probes, not an exhaustive compressed-offset check. The full scan
streams lines (memory is bounded by a FASTA line and metadata, not the total genome size). Marker
queries use the same index checks without requiring a Picard dictionary.
Canonical dictionary naming is shared here instead of duplicated in shell and Python.

### prepare_y_refs.sh
Stages hs1 and hg38 reference FASTAs into `./resources/ref/` and ensures
required sidecars exist (`.fai`, `.dict`). By default it symlinks source
FASTA files; use `--copy` to copy. Use `--primary-only` to build staged
FASTA files containing canonical contigs only (`chr1-22,chrX,chrY,chrM` by default).
Staging accepts relative source paths, preserves BGZF `.gzi` and canonical dictionary sidecars,
and does not replace an existing target except to reuse an unchanged source symlink. Use a fresh
output directory for a new copy/subset; orphaned sidecars are rejected. Subsets get fresh indexes
for their new byte layout rather than copied source offsets. Failed extraction does not leave a
partial target FASTA.
The default GRCh38 staging name is now `hg38.fa`; use `--hg38-name` if another label is needed.
`Y_RESOURCES_DIR` changes the default staging root to its `ref/` subdirectory; `--out-dir` wins.

### run_y_dual_liftover_experiment.sh
Two-branch experiment driver for modern Y analysis:
- branch A: liftover markers `hg38 -> hs1`, then call against native hs1 sample
- branch B: liftover sample `hs1 -> hg38`, then call against hg38 markers
- branch C (if GFF3 markers are available): rank hg38-lifted sample with hg38 marker GFF3 (`yfull_node` labels) for terminal-style naming
- branch D: path ranking with auto top-level clade (default) or manual `--clade-prefix`
- branch E: condensed `top_paths.tsv` summary (top positive-score rows per dataset)

The optional sample gVCF is queried on the native hs1 branch only, not lifted. GFF3 interpretation
is an additional hg38 branch, so different panels can yield different evidence. Rerunning only one
ranker leaves other reports stale; save matching run/script/resource revisions when comparing.
Each extraction uses its matching hs1/hg38 reference through the common resolver; pysam is required.
Optional `--min-gq/--min-dp` apply to all branches' marker states (default: unset).

Default output path is under the caller's working directory: `./experiments/...`.
Run it from a private working directory for sensitive sample data.
Supports `--java-opts` (or env `PICARD_JAVA_OPTS`) for LiftoverVcf heap sizing.
Uses `RECOVER_SWAPPED_REF_ALT=true` to retain loci where reference/alternate
alleles are swapped between assemblies.
Accepts `--markers-gff3-hg38` to force a specific GFF3 marker file.
References can be passed explicitly from the index build; prefer the compatible base assembly to
large appended mapping references. The script validates reference sidecars and indexed access, but
does not establish whole-sequence identity by default; sidecar/probe checks and optional full MD5 validation
are described above.
If refs are omitted, shared discovery searches:
- `--ref-root` or `REF_ROOT`/`BIO_TOOLS_INDEX_DIR`/`INDEX_DIR`
- `$Y_RESOURCES_DIR/ref` when configured
- `./resources/ref`
- `mapping/index` under the repo root
Base analysis references are preferred across these roots, even when extended fallbacks are allowed.
By default the script rejects references that look like extended H3630/decoy sets;
use `--allow-extended-ref` only when you explicitly want that behavior.

Hard-won operational notes:
- For routine Y analysis, prefer non-extended refs (`chm13v2.0_maskedY_rCRS.fa(.gz)` and
  `GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna`) to avoid decoy-driven confusion.
- Picard liftover is strict about reference sidecars; shared validation now checks their contents.
- For `ref.fa.gz`, the canonical dictionary is `ref.dict`; the helper aliases the previously
  miscomputed `ref.fa.dict`. Do not infer validity from the presence of either filename.
Canonical usage:
- Male sample QC + assignment cross-check: use this script.
- It is the maintained modern Y entrypoint.

### liftover_to_hg38_batch.sh
Batch liftover helper for ordinary hs1/CHM13 sample VCF files to GRCh38.
Designed for generating GRCh38 callsets for downstream tools (Exomiser/PharmCAT).
The shell entry point delegates to `liftover_to_hg38_batch.py`; there is one batch implementation.
Full gVCF conversion is **unsupported**: a streaming content check rejects `<NON_REF>`, `<*>`,
and nonvariant `ALT=.` records with `END`, regardless of filename or gzip suffix. Unused header
declarations alone are allowed; literal spanning-deletion `*` is not the `<*>` confidence allele.
See the [Picard gVCF limitation](LIFTOVER.md#picard-gvcf-limitation).
Optional `--primary-only` emits a canonical-contig subset
(`chr1-22,chrX,chrY,chrM`) per lifted sample.
Supports `--java-opts` (or env `PICARD_JAVA_OPTS`) for LiftoverVcf heap sizing.
Uses `RECOVER_SWAPPED_REF_ALT=true` to retain loci where reference/alternate
alleles are swapped between assemblies.
Canonical usage:
- For non-Y pipelines (including female exome), use this script for hs1 -> hg38 conversion.
- `--ref-hg38` can be omitted if your index root is discoverable via `--ref-root` or `BIO_TOOLS_INDEX_DIR`.
  `REF_HG38` and trusted local config are also supported. Extended reference filenames now require
  `--allow-extended-ref`, matching the dual-Y driver; their presence does not authorize auto-selection.
- Supply a new or empty dedicated `--out-dir`; populated directories and destination symlinks
  are refused. All inputs are scanned before reference-sidecar creation or output/report writes.
  Duplicate inputs (including hardlinks), colliding stems and malformed/truncated streams fail early.
- Paths containing spaces and dotted stems are preserved. An ordinary variant-only file with a
  `.gvcf` suffix is treated as VCF, not as a validated gVCF conversion; same-stem inputs collide.
- Accepted/rejected VCFs, indexes, optional primary subsets, `liftover_report.tsv`, `manifest.json`
  and tool logs are staged in a hidden sibling directory. The batch is published with one directory
  rename only after decoding, indexed reads, sample-header checks and record conservation pass.
  `--primary-only` still filters after conversion, not before the source scan.
- Failed/interrupted stages are retained with logs and a failure manifest, without partial final
  products. No automatic resume or overwrite is provided. Run locks prevent cooperating runs from
  sharing a destination; inspect stale locks after a crash rather than blindly deleting them.
- Requires Python/pysam plus Java/Picard, samtools, bcftools and tabix. Batch lifecycle regressions:
  `python3 annotate/test_liftover_batch.py` (also in `sh annotate/run_tests.sh`). These use a stand-in
  Picard process and real indexing/subsetting tools, not production sample conversion.

See [batch lifecycle and limits](LIFTOVER.md#what-the-existing-drivers-actually-do) for recovery and
the difference between atomic publication, reference metadata changes and power-loss durability.

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
- `<prefix>.matrix.tsv`: original status/GT/DP/GQ/HG columns plus per-dataset coordinates,
  marker alleles, source and reason; ID-mode rows therefore retain both assemblies' positions.
- `<prefix>.summary.tsv`: per-dataset status counts and the legacy global-union agreement/flip
  table; optional expected-panel totals precede that final table.
- `<prefix>.pairwise.tsv`: explicit global/pair unions, shared keys, jointly callable counts,
  callable agreement percentage, distinct callable locus-pair counts and allele-check diagnostics.
- `<prefix>.expected.tsv`: only with `--expected-states`, one row per expected ID and dataset,
  including missing panel markers, filtered calls and ancestral/derived conflicts.

Use `--key-mode locusid` only on a common assembly; `--key-mode id` needs verified marker identities
across assemblies. Duplicate dataset names, malformed rows and duplicate comparison keys are
rejected before reports are opened. `--duplicate-policy collapse-identical` permits only exact
duplicate rows, reports the collapse and never chooses between conflicting calls. Conflicting
alias calls for the same locus and ancestral/derived contrast are also rejected.

Allele checks use the same ancestral/derived rule as the resolver (`marker_schema.py`). REF/ALT
swaps preserving biological alleles are allowed. Complemented **single-base** definitions require
`--key-mode id --allow-reverse-complement`; this is not proof of cross-assembly coordinate identity.
Incompatible definitions fail even if their reported statuses agree. Callable rows with unknown
allele definitions are rejected; unresolved rows remain visible as unverified evidence.

`callable_both` includes only derived/ancestral states on both sides. Filtered, missing, ambiguous,
unknown and other states do not enter that denominator; zero callable pairs give `NA`, not 100%.
Counts are per comparison key and may include aliases; `callable_unique_locus_pairs` exposes the
number of distinct coordinate pairs, not independent molecules. All-row agreement still includes
matching no-calls and missing/missing pairs and must not be presented as accuracy.

For expected-path QC use a user-supplied TSV with `id` and `status` columns, where status is
`derived` or `ancestral`. The full supplied panel is the denominator, even when an ID is absent.
Ambiguous multi-locus IDs and expected IDs excluded by other filters are rejected. The tool does
not infer a path, decide equivalence, or validate the expected panel as ground truth. Expected
totals count supplied IDs, not independent loci: equivalent aliases can represent the same allele. Keep modern
expected panels and every per-sample report outside the repository.

```bash
python3 "$REPO_ROOT/annotate/compare_marker_status_sets.py" \
  --input "native=$ANALYSIS_DIR/native.marker_status.tsv" \
  --input "lifted=$ANALYSIS_DIR/lifted.marker_status.tsv" \
  --key-mode id --expected-states "$ANALYSIS_DIR/expected-path.tsv" \
  --out-prefix "$ANALYSIS_DIR/comparison"
```

Existing reports are refused unless `--force` is explicit; input files can never be overwritten.
An old expected-state report requires the matching option or a fresh prefix to prevent stale
partial reruns. Validation precedes report creation, but later I/O failure is not transactional.
AD and definitive liftover-loss provenance are not available in the current input schema; see the
[comparison guide](LIFTOVER.md#comparing-native-hs1-lifted-hg38-and-direct-hg38) for remaining limits.

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
