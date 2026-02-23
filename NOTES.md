# Notes

**Dataset**
Cell Genomics 2023 short article (Sep 13, 2023): "High-coverage genome of the Tyrolean Iceman reveals unusually high Anatolian farmer ancestry." This note set refers to the Otzi 2023 resequencing BAMs.

**Path Conventions**
- `<raw_data_dir>`: location of original input BAM files.
- `<analysis_dir>`: location of per-run working/output files.
- `<legacy_data_dir>`: location of legacy comparison VCF datasets.
- `<output_dir>`: location of final variant-calling outputs.

**Findings**
- EAGER/AdapterRemoval BAMs can be merged and mate-stripped. MarkDuplicates will treat singletons as single-end and can over-mark duplicates in repetitive regions.
- Pseudo-pairing merged reads enables BWA MEM insert-size modeling. BWA MEM typically reports only FR orientation for Illumina/BGI libraries because other orientations are too rare to model.
- For remap experiments, the safest baseline is the authors' final merged and deduped BAM. Re-dedup can bias results when mates are missing.
- `samtools fixmate` does not create pairs; it only fills mate fields for reads that are already correctly paired.
- Paper QC: contamination estimated with ANGSD from haploid X heterozygosity. Reported 0.5% ± 0.06% contamination vs 7.5% ± 0.25% in the earlier genome.
- `bwa-postalt.js` empty-line guard (`file.readline(buf) > 0`) is upstream since 2018-04-07 (bwa commit `27dd1da`); keep upstream version.
- `samtools fastq` collapses multiple READ_OTHER records with the same QNAME down to a single record. If F/R reads lack READ1/READ2 flags, the R read can disappear unless you set flags first.
- RevertSam(Spark) with `SANITIZE=true` drops/normalizes problematic mate records. This breaks F/R flag semantics and can cause `samtools fastq -1/-2` to lose mates.

**Troubleshooting (Regression): MarkDuplicatesSpark Pairing Crash**
Error:
`GATKException: Detected multiple mark duplicate records objects corresponding to read with name '...'`

Symptoms:
- Same QNAME appears twice with flags `0`/`16` in BWA output (unpaired).
- MarkDuplicatesSpark aborts on those records.

Checks:
- `samtools view -c -f 1 bwamem.bam` (paired count)
- `samtools view -c -F 1 bwamem.bam` (single-end count)
- `samtools view -e 'qname=="<QNAME>"' bwamem.bam` (or `grep "^<QNAME>"`)
- `samtools fastq -n ubam.bam | head -8` (verify interleaved pairs are adjacent)

Potential causes (use as a checklist):
- Fastq stream de-interleaved (e.g., fastp run with `--interleaved_in` but without `--interleaved_out`).
- Pair adjacency broken by filtering or reordering steps between fastq and BWA.
- Paired input not passed to BWA (interleaved fastq without `-p`, or paired fastqs not used).

Fix pattern:
- Ensure the stream remains interleaved to BWA (`bwa mem -p`), or emit paired fastqs and run `bwa mem ref r1 r2`.
- After pipeline fixes land, keep this section as a regression checklist for unexpected pairing loss.

**BWA Pairing Behavior (Source Notes)**
- `bwa mem -p` uses `bseq_classify()` to split reads into SE vs PE by checking *adjacent* read names for equality.
- Read names are normalized by `trim_readno()` (suffix `/1` or `/2` is removed before comparison).
- No BAM tags or flags are used at this stage (FASTQ has no flags). Pairing is purely name + adjacency.
- Consequence: any step that de-interleaves or reorders the stream (e.g., fastp without `--interleaved_out`) turns true pairs into single-end reads.
- Using paired FASTQs (`-1/-2`) bypasses this classification; singletons may be dropped depending on `samtools fastq` flags.
- With singletons present, `bwa mem -p` can still split some true pairs across chunk boundaries (it only preserves an even *count* of reads per chunk, not pair boundaries). Those pairs get treated as single-end (flags 0/16). Mitigations: use paired FASTQs (`-1/-2` with `-s` for singletons), or set a larger `-K` chunk size to reduce splits.

**Pipeline Stages Where Pairing Can Break**
- `eager_repair_bam.py`: can create or drop pairs if run in the wrong mode.
- uBAM creation: should preserve pairing flags in BAM, but FASTQ export discards flags.
- `samtools fastq`: interleaving/paired output options control adjacency; singletons can be dropped.
- `fastp`: must preserve interleaving if used in a stream.
- RevertSam(Spark) with `SANITIZE=true`: can remove mates or clear pairing semantics, leading to READ_OTHER collapse in fastq.

**samtools fastq READ_OTHER Gotcha**
- `samtools fastq -1/-2` only routes reads with READ1/READ2 flags.
- Reads with neither (READ_OTHER) go to `-0`. If both mates are READ_OTHER, only one is output (QNAME collapse).
- Fix: set READ1/READ2 flags first (`eager_repair_bam.py --pair` on name-sorted input), or split by XQ tag.

**DeDup Behavior (Source Notes)**
- Expects coordinate-sorted input; it does not sort. Unsorted input can miss duplicates.
- Output header is set to unsorted; `-u/--unsorted` flag is effectively broken (always forces unsorted header).
- DeDup removes duplicates (does not mark with 0x400). Existing duplicate flags are preserved on surviving reads.
- Requires `M_/F_/R_` prefixes unless `--merged` is used.
- DeDup ignores read-group (RG). Original protocol: DeDup each library, then merge and DeDup again (likely to reduce compute and catch cross-library duplicates).
- Per-library DeDup also supports per-library QC reporting (duplication rates by library); final DeDup captures cross-library duplicates after merge.

**eager_repair_bam Usage**
Tags:
- `XQ:Z:M` original merged read
- `XQ:Z:D` synthetic duplicate mate
- `XQ:Z:F/R/U` original prefix or unknown

Examples (fast BAM pipes, no SAM I/O):
```bash
# Strip prefixes + duplicate merged reads, then name-sort
samtools view -u input.bam \
  | ./util/eager_repair_bam.py --strip-prefix --duplicate-merged --dup-orientation flip -o - --output-uncompressed \
  | samtools sort -n -@8 -m 2G -o stripped.qn.bam -

# Pair reads in name-sorted BAM (uses XQ tag)
./util/eager_repair_bam.py --pair stripped.qn.bam -o paired.bam --output-uncompressed

# Optional mate-field fill (needed for samtools markdup)
samtools fixmate -@8 -m paired.bam paired.fixmate.bam
```

Useful filters:
```bash
# Synthetic duplicates only
samtools view -c -e 'XQ=="D" && (flag&0x400)' file.bam

# Non-synthetic duplicates
samtools view -c -e 'XQ!="D" && (flag&0x400)' file.bam
```

**Current Run Results (Otzi D2049 combined)**
- Baseline published BAM prefix counts:
  - `M=138980507`
  - `F=5545558`
  - `R=3582337`
- Reconstructed checkpoint (`pair.bam`) flag sanity:
  - total `287088909`
  - read1/read2 `144526065 / 142562844`
  - singletons `5182509`
- Primary-assembly comparison:
  - pre-DeDup (`pair.prim.bam`): `145064186` total, `M/F/R=136687692/5138191/3238303`
  - post-DeDup (`pair.prim_rmdup.bam`): `138328921` total, `M/F/R=131715629/3743710/2869582`
- DeDup log:
  - reverse removed `368721`
  - forward removed `1394481`
  - merged removed `4972063`
  - total removed `6735265`
  - reported duplication rate `0.05` (rounded), exact `4.64%`
- Loss rate by class:
  - `M: -3.64%`
  - `F: -27.14%`
  - `R: -11.39%`

**Run Interpretation**
- DeDup did not behave like insert-size/template-aware PE deduplication.
- Loss was concentrated in unmerged classes (`F/R`), especially `F`.
- Keep this DeDup result as one analysis branch, and document singleton-driven over-pruning risk in conclusions.

**Final Merged DeDup Pass (D2049-D2052)**
- Pre-final-DeDup (`...sort.bam`) total: `668715804`
- Post-final-DeDup (`...sort_rmdup.bam`) total: `656719414`
- Removed in final pass: `11996390` (`1.79%`)
- Mapping rate remained stable at `99.98%` before/after.

**DeepVariant Run Notes (Ancient DNA Branch)**
- DeepVariant WGS run is in progress on primary-assembly final BAM.
- `call_variants` throughput stabilized around `~2.66-2.70 sec per 100` examples on CPU-only execution.
- `--disable_small_model=false` is intentional for default performance.
- One full run completed successfully (`iceman.vcf` / `iceman.gvcf`).
- A second side run (with `--regions`) was interrupted by power loss; usable resume state depended on whether intermediate files were persisted outside container `/tmp`.
- Observed long `postprocess_variants` tail was active compute (protobuf `_message.so` hot path), not idle hang.

**DeepVariant Reference Compatibility Gotcha**
- Error `Reference contigs ... IS MISSING` with `chr1..chrY` usually indicates BAM/reference mismatch, not "too many contigs".
- Confirm by comparing BAM `@SQ` lengths to FASTA `.fai`.
- Example observed mismatch: CHM13-aligned BAM lengths vs GRCh38 reference lengths.

**Y Haplogroup Re-check (Published vs SOLiD vs Illumina)**
- Literature baseline used for comparison:
  - Nature Communications 2012 (`ncomms1701`): Iceman assigned to haplogroup `G`, subgroup `G2a`, including `M201` and `L91` context.
  - Nature Communications 2025 (`s41467-025-61601-8`): reports lineage `G2a2a1a2a1a1b (G-Z6208*)`.
- Legacy SOLiD-derived branch (`<legacy_data_dir>/chrY_raw_Iceman_tst_hg38.vcf.gz`) remains noisy and conflict-heavy with current marker panel:
  - broad cross-clade signal and no consistent deep G-branch support.
  - key markers often unresolved/no-call under stricter filtering.
- Illumina reanalysis branch (DeepVariant output `<output_dir>/iceman.vcf`) shows coherent G-branch signal:
  - `M201` derived (`GT=1/1`, `GQ=42`, `DP=10`)
  - `L91` derived (`GT=1/1`, `GQ=23`, `DP=7`)
  - `L166` derived (`GT=1/1`, `GQ=14`, `DP=3`)
  - `Z6208` derived (`GT=1/1`, `GQ=39`, `DP=9`)
- Current interpretation:
  - Illumina branch is consistent with the published `G2a` direction and supports a downstream `Z6208*` placement.
  - SOLiD branch is not sufficient alone for robust terminal placement with the current marker set.
- DeepVariant small-model on/off comparison on Iceman:
  - all-sites shared GT differences are non-zero (`~0.76%`), but PASS-only shared GT differences are low (`~0.017%`).
  - haplogroup-defining markers (`M201`, `L91`, `L166`, `Z6208`) stayed `PASS` + `GT=1/1` in both runs.
  - `Z6208` confidence dropped in no-small run (`GQ/QUAL` lower), but call direction did not change.

**Y Marker Analysis Tooling Added**
- `annotate/y_haplo_from_vcf.sh`
  - marker-merge + derived-call extraction workflow (YBrowse hg38 marker VCF).
  - caller-aware filter modes: `any|pass|pass-or-dot|deepvariant`.
  - supports DeepVariant outputs with nontrivial FILTER semantics and gzip inputs without `.gz` suffix.
- `annotate/y_haplo_from_markers.py`
  - marker-state extraction for both VCF and gVCF inputs (`derived|ancestral|nocall|ambiguous` per marker).
  - intended to preserve no-call/ancestral evidence, not only derived rows.
- `annotate/y_clade_consistency.py`
  - compares candidate clades across datasets using support/conflict counts.
  - intended for published/SOLiD/Illumina side-by-side scoring instead of raw "deepest label string" heuristics.
- `annotate/y_path_rank.py`
  - root-to-tip style clade ranking from marker-status rows with tunable scoring.
  - includes optional down-weighting of potential deamination transitions (`C>T`, `G>A`) in derived evidence.

**Method Attribution**
- Path-oriented haplogroup scoring and damage-aware marker handling are aligned with practices described in:
  - https://doi.org/10.1101/2024.03.13.584607
- Parsimonious path interpretation with ancestral-state evidence is conceptually aligned with:
  - PathPhynder (as cited by that work).

**YBrowse Marker Source Notes (2026-02-23 refresh)**
- `snps_hg38.vcf.gz` and `snps_hg38.gff3` are not row-equivalent exports.
- Observed differences are mainly due to representation:
  - GFF includes indels and more duplicate loci rows.
  - VCF export is SNP-oriented and has fewer duplicate-locus rows.
- Practical workflow decision:
  - Use marker VCF as the canonical input for current haplogroup scripts.
  - Sanitize marker VCF alleles before liftover/strict parsing (`ALT!=REF`, no duplicate ALT).
