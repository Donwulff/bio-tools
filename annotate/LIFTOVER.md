# Liftover and Cross-Assembly QC

This guide records the reusable conclusions from the hs1/CHM13v2-with-Y and GRCh38 experiments.
Modern sample files, identifiers and results belong in an external analysis directory. For public
Iceman results and overall status, start at [the situation review](../SITUATION.md).

## Separate Conversion From Y Validation

`liftover_to_hg38_batch.sh` converts sample callsets regardless of sex or assay. Y analysis is an
optional validation on informative male samples; it is not a prerequisite for exome conversion.
`run_y_dual_liftover_experiment.sh` tests both coordinate directions: hg38 markers lifted to hs1
against native calls, and sample variants lifted to hg38 against hg38 markers. Its optional sample
gVCF is read on the native hs1 branch; that driver does **not** lift the sample gVCF to hg38.

The shared implementation is currently `lib_liftover.sh` for tool/dictionary helpers, plus
`y_haplo_from_markers.py` for marker states, `y_path_rank.py` for descriptive ranking, and
`ytree.py` for an explicit topology where one is provided. Reference discovery and Picard invocation
still occur in multiple drivers. Compatibility warnings on the older scripts do not remove that
duplication. A future consolidation should reuse these operations rather than add another caller.

## References and Resources

- Keep the full competitive reference used for **mapping** distinct from the smaller target FASTA
  used for **liftover**. A chrY-only marker query does not justify chrY-only read alignment; the latter
  loses competing X/autosomal loci. See the [recovery-versus-specificity correction](../NOTES.md#the-chry-only-finding-measured-recovery-not-contamination-2026-07-27).
- Prefer the already-downloaded base `chm13v2.0_maskedY_rCRS.fa(.gz)` and an appropriate GRCh38
  analysis FASTA for liftover. The index builder also downloads
  `GCA_000001405.15_GRCh38_full_analysis_set_masked.fna.gz`; the no-alt version is another explicit
  choice. Do not silently substitute appended HLA/IMGT, oral or other mapping additions.
- Check the chain's assembly, direction, sequence names and target sequence content. A custom
  reference does not by itself require new chains if the relevant named sequences are unchanged;
  altered sequence, masking or mitochondrial choice requires locus-specific compatibility checks.
- Use explicit reference/marker/chain paths when reproducing a run. `BIO_TOOLS_CFG`,
  `BIO_TOOLS_INDEX_DIR`/`INDEX_DIR` and `Y_RESOURCES_DIR` are conventions with uneven driver support,
  not a fully implemented common resource registry. Running from another directory changes some
  fallback paths and can trigger a second download.
- Save source URLs, retrieval/header dates, hashes, command lines and tool versions alongside a run.
  A copyright year is not a release date, and an old local fallback is not a successful refresh.
- YBrowse VCF and GFF3 are different exports. GFF3 carries ancestral/derived alleles and YFull
  labels missing or coarse in some VCF rows; unsorted rows and aliases need explicit handling.
  The current GFF3 resolver selects single-position SNP/point rows, so fetching GFF3 does not imply
  general indel support. Sanitize invalid VCF ALT definitions before strict Picard parsing.

Reference preparation needs a readable nonempty FASTA, usable `.fai`, and matching nonempty
dictionary containing `@SQ`. Compressed FASTA suffixes do not establish BGZF compatibility. For
random access, verify BGZF and its `.gzi` or stage an uncompressed FASTA. A `.gzi` is not a VCF
tabix index. Do not reuse byte-offset indexes after copying, decompressing or rewriting sequence
layout; regenerate indexes for the exact staged bytes.

For `x.fa.gz`, the canonical dictionary is `x.dict`. The shared helper also maintains the earlier
`x.fa.dict` alias. The original failure came from the script's suffix calculation; there was no
demonstrated need to explain it as inconsistent Picard behavior. Current helpers check existence,
not dictionary contents or agreement with the FASTA. An empty/stale sidecar can therefore still
pass their preflight. `prepare_y_refs.sh` links/copies by default; `--primary-only` is required to
extract canonical contigs. Its output filename alone does not prove a subset was made.

## Picard gVCF Limitation

The earlier session explanation that millions of `IndelStraddlesMultipleIntevals` rejects were
merely expected chain-boundary losses was incomplete. **Picard 2.27.1 cannot be assumed to preserve
gVCF reference-confidence blocks.** Its source compares REF length with lifted interval length.
A block with a one-base REF and a multi-base `END` interval can fail that comparison even when
the interval maps without a gap. The filter's name therefore does not prove an indel or chain
discontinuity. Picard buffers accepted records for sorting while writing rejects during the scan,
which also explains delayed accepted output. See the pinned
[LiftoverVcf source](https://github.com/broadinstitute/picard/blob/2.27.1/src/main/java/picard/vcf/LiftoverVcf.java).

Keep native gVCFs for reference/no-call evidence. Treat existing lifted gVCFs as experimental
partial outputs; neither successful exit nor some accepted variants validates their block semantics.
Use a separately validated block-aware transformation, direct target-reference calling, or targeted
evaluation of lifted marker sites on the native assembly when reference-confidence evidence is
required. Lowering minimum match does not repair block semantics. No such replacement conversion
was validated in this review.

For ordinary VCF rejects, inspect FILTER, variant type, contig and source-to-target representation
separately. `NoTarget`, allele mismatch and interval-length rejection are different failure classes.
Reject coordinates remain in the source assembly. REF/ALT recovery must preserve the biological
ancestral/derived interpretation; auxiliary allele annotations need checking after allele swaps.

## What the Existing Drivers Actually Do

The batch script now distinguishes `.vcf` and `.gvcf` outputs and refuses to overwrite those output
files. This repairs the original collision when both inputs had the same sample stem. The check
does not establish transaction safety: it writes the run report before all per-file checks, and a
failed run can leave products or sidecars. Use a fresh output directory and retain the full log.
`--primary-only` selects target canonical contigs **after** liftover; it does not restrict the
input scan or remove source-decoy processing beforehand.

Example template for ordinary VCF conversion (paths are user configuration):

```bash
REPO_ROOT=/path/to/bio-tools
ANALYSIS_DIR=/path/to/external-analysis
REF_HG38=/path/to/compatible-hg38.fa
CHAIN_HS1_TO_HG38=/path/to/hs1ToHg38.over.chain.gz

"$REPO_ROOT/annotate/liftover_to_hg38_batch.sh" \
  --chain-hs1-to-hg38 "$CHAIN_HS1_TO_HG38" \
  --ref-hg38 "$REF_HG38" \
  --picard-jar /path/to/picard.jar \
  --java-opts "-Xms1g -Xmx8g" \
  --out-dir "$ANALYSIS_DIR/liftover-run" \
  --primary-only \
  --sample "$ANALYSIS_DIR/input.vcf.gz"
```

This reflects current CLI options; it is not a claim that those paths or that heap fit every input.
Reference loading and VCF sorting are separate memory costs. The Picard calls in these drivers use
the detected jar, even if GATK is available for dictionary creation; record the actual backend.

## Comparing Native hs1, Lifted hg38 and Direct hg38

First match marker definitions, alleles, site filters and quality thresholds. Keep failed marker
liftover, missing sample records, explicit no-call, confident reference/ancestral calls and QC-rejected
calls separate. A variant-only VCF's missing row is not evidence for the reference allele. Native
gVCF helps, but its block-start REF alone does not give the base at an interior marker: resolving
that state needs the sample assembly's reference sequence. The resolver now takes `--reference`
and performs indexed lookup at the marker coordinate, checking the record's anchor REF against
the FASTA. A `0/0` block can be **derived** relative to a catalogue ancestral allele: reference
genotype is not synonymous with ancestral. Both `<*>` and `<NON_REF>` blocks are recognized.

For this lookup, install pysam and prepare the native FASTA's nonempty `.fai`; BGZF also needs
`.gzi`. Plain gzip is not a substitute. The resolver does not create reference indexes or load a
whole chromosome into memory. A missing reference leaves interior markers `nocall` with reason
`reference_required`; incompatible/masked bases remain no-calls with a diagnostic reason. Overlapping
reference blocks are rejected rather than arbitrarily choosing an interval. Block markers are
currently limited to single-base allele definitions, not indels or general symbolic variants.

`--min-gq/--min-dp` now filter before **all** state outputs and summaries, not just derived reports.
Failures are `filtered`, with raw GT/GQ and effective DP retained. Excluded FILTER records are kept
as filtered evidence, not discarded into "missing" or replaced by a covering reference block.
The appended 14th column records the reason; the original 13 columns retain their meaning/order.
The ranker treats filtered rows as no-calls (default weight zero); comparisons retain their status.
For blocks, effective depth prefers `MIN_DP` when present and falls back to DP only when MIN_DP
is absent. This is a conservative block bound, not an exact per-marker allele depth. Zero effective
block depth cannot supply a reference/ancestral call, even with thresholds unset. The block and
MIN_DP conventions follow the [GATK gVCF description](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format).

Example for a **native** gVCF and a panel on that same assembly:

```bash
python3 "$REPO_ROOT/annotate/y_haplo_from_markers.py" \
  --input "$ANALYSIS_DIR/input.gvcf.gz" \
  --markers "$MARKERS_NATIVE" --reference "$REF_NATIVE" \
  --site-filter-mode deepvariant --min-gq 20 --min-dp 3 \
  --out-prefix "$ANALYSIS_DIR/y-qc"
```

Thresholds here are illustrative, not registered research changes. Synthetic regressions cover
lookup and filtering; existing sample outputs have not been regenerated or validated by this fix.
The dual driver passes the appropriate native/lifted reference automatically; the Iceman comparison
and compatibility modern drivers accept `--reference`. Preserve old outputs and use a fresh run
directory, especially because the Iceman driver's existing default thresholds now affect ranking.
This fixes native marker interpretation, **not Picard gVCF block liftover**.

`compare_marker_status_sets.py --key-mode locusid` is for a common coordinate system. Across
assemblies use verified marker identities with `--key-mode id`, or a validated coordinate/allele
crosswalk. The comparator now rejects conflicting/duplicate keys, contradictory aliases and
incompatible biological allele definitions. Exact duplicate rows can be collapsed only by an
explicit policy. Complemented SNP definitions need explicit ID-mode permission and still require
an independently verified identity/coordinate crosswalk; this is not general indel normalization.

Use the new `.pairwise.tsv` callable-on-both denominator alongside the retained legacy all-row
summary. Only derived/ancestral pairs are callable; zero callable evidence yields `NA`. Matching
no-calls or missing/missing rows can still dominate all-row agreement. The legacy global-union
counts also include keys seen only in a third dataset; global and per-pair unions are now explicit.
Per-key counts may include aliases; distinct coordinate-pair counts expose that duplication without
claiming to count independent molecules. Neither agreement statistic establishes accuracy.

Optional `--expected-states` takes a trusted ID/state panel and reports the entire panel per dataset,
including IDs missing from the input. Filtered/ambiguous/no-call states remain unresolved, not
ancestral conflicts. This is an explicit expected-panel comparison, not automatic tree/path
construction or proof of the expected states. Missing marker rows alone do not prove liftover
rejection: unqueried panel entries and rejected coordinates still need separate provenance.
The input schema has no allele-specific AD, so the comparator does not invent supporting-read
counts. See [the comparator options](README.md#compare_marker_status_setspy).

Within a trusted expected Y path, review ancestral/derived conflicts, equivalent-marker support,
DP, allele-specific AD, GQ and BAM-level MAPQ/strand evidence. Count unique loci and molecules rather
than every alias independently. Current `derived_dp_sum/mean/max` describe depth at derived marker
rows; they neither select the derived component of AD nor deduplicate molecules shared across loci.
The default deamination weight is a heuristic inherited from aDNA and is not a probability or an
appropriate automatic penalty for every modern transition.

The three callsets share underlying reads and are not independent truth sets. A known terminal
assignment supplies a limited expected-path control, not genome-wide ground truth. Neither direct
GRCh38 nor lifted CHM13 is automatically the winner; inspect discordant loci and callable regions,
and use independent truth data when making accuracy claims. External annotation tools' supported
assembly/contig sets must be checked for the particular tool and release.
