# Notes

**Dataset**
Cell Genomics 2023 short article (Sep 13, 2023): "High-coverage genome of the Tyrolean Iceman reveals unusually high Anatolian farmer ancestry." This note set refers to the Otzi 2023 resequencing BAMs.

**Path Conventions**
- `<raw_data_dir>`: location of original input BAM files.
- `<analysis_dir>`: location of per-run working/output files.
- `<legacy_data_dir>`: location of legacy comparison VCF datasets.
- `<output_dir>`: location of final variant-calling outputs.

**Data Policy**
- Ancient/public sample results are in-scope for repository docs and commits.
- Modern/private sample IDs and per-sample outputs are out-of-scope for committed docs; keep them in external analysis directories and use placeholders in tracked notes.
- Naming gotcha: `.gitignore` has a broad `*private*` guard for personal analysis outputs. In Y phylogenetics "private variant/branch" is standard terminology, so a legitimately committable file can be silently ignored. `results/iceman_y_private_branch_candidates.tsv` was invisible to git for this reason; renamed to `..._novel_branch_candidates.tsv`. Prefer "novel" in filenames and leave the guard alone.

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

**EAGER aDNA BAMs → Modern GRCh38 Pipelines (Iceman / Otzi re-analysis findings)**
Published EAGER-processed BAMs are common in aDNA. They are typically:
- Merged reads (AdapterRemoval) with `M_ / F_ / R_` prefixes.
- Mate information stripped for unmapped or filtered reads.
- Already deduplicated against an older reference.

Feeding them to modern pipelines (revert + BWA + MarkDuplicates + variant calling on GRCh38/hs1) is mostly straightforward using the repair tool, **but**:
- Standard Picard/GATK `MarkDuplicates` (even with pseudo-paired reads via `eager_repair_bam.py`) falls back to single-end logic on stripped data and over-marks in repetitive regions.
- DeDup is required for correct aDNA-aware deduplication on merged reads and to produce comparable "passing reads" counts to the original paper.
- The goal of the re-analysis was to enable fair comparison on updated references while preserving the ability to match published metrics.

Recommended reusable flow (see also `util/eager_repair_bam.py` and `mapping/eager_prep.sh`):
1. Start from the published EAGER BAM.
2. Name-sort + strip prefixes + duplicate merged reads (for PE-tool compatibility):
   ```sh
   samtools sort -n input.bam -o qn.bam
   ./util/eager_repair_bam.py --strip-prefix --duplicate-merged qn.bam -o repaired.bam
   ```
3. `samtools fixmate` + coordinate sort.
4. For DeDup branch: re-apply prefixes (or keep the tagged version) and run DeDup on coordinate-sorted input.
5. Keep both the "standard pipeline" output and the DeDup branch for comparison.

**Building DeDup (third-party Gradle/Java tool)**
The copy under `DeDup/` here is only for reference (the entire tree is ignored via `.gitignore` in this repo).

To build:
- Clone the maintained fork: `git clone https://github.com/Donwulff/DeDup.git`
- `cd DeDup && ./gradlew jar`
- The fat jar lands in `build/libs/DeDup-0.12.9.jar`

Requires a JDK 8+ (with javac). On some Linux setups ensure the `-jdk` package (not just JRE) is installed for the desired Java version.

Other options: bioconda (`conda install -c bioconda dedup`) or download a pre-built jar from the fork's releases.

See the fork's README for usage. The wrapper + build.gradle fixes live in the fork (not vendored long-term here).

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
  - Illumina branch is consistent with the published `G2a` direction and confirms `Z6208` derived.
  - Terminal placement was later resolved at read level as `G-L166*`, **not** a node downstream of `Z6208`
    — see "Terminal Y Placement from Read-Level Evidence" below. `Z6208` is an L166-level SNP.
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
  - For modern terminal-label comparisons against YFull-style nodes, marker GFF3 can be used in `y_haplo_from_markers.py` because it includes `yfull_node` labels.
  - YBrowse GFF3 is not position-sorted; marker evaluation must not assume monotonic marker order.

**Liftover Ops Gotcha (2026-02-24)**
- A frequent failure mode was not the chain or markers, but reference sidecar naming for compressed FASTA.
- For `.fa.gz` references, `CreateSequenceDictionary` may produce `ref.fa.dict` while `LiftoverVcf` path checks may look for `ref.dict` (tool/path dependent).
- Mitigation now used in scripts:
  - create/use the canonical dictionary path based on FASTA suffix,
  - add dictionary alias symlink so both names resolve,
  - preflight `.fai` + `.dict` before liftover starts.
  - reuse shared liftover helper module (`annotate/lib_liftover.sh`) so this behavior is consistent across drivers.

**Y Path Ranking Readability/Label Gotcha (2026-02-24)**
- Some marker tables carry coarse clade labels (e.g. `I1`) even when marker IDs are specific (e.g. `Y47125`).
- Without ID-aware labeling, path outputs can show specific SNP labels as `score=0` even when the marker is derived in status tables.
- Current ranking adds synthetic `<top-clade>-<marker_id>` labels for these rows (example: `I-Y47125`), and dual-liftover runs now emit compact `top_paths.tsv` summaries.

**Terminal Y Placement from Read-Level Evidence (2026-07-25)**
Direct pileup on the final aDNA BAM (`...pair.prim_rmdup.sort_rmdup.coord.bam`), not caller output,
resolves the terminal node. Filters: `-q 25 -Q 20 --no-BAQ`, hg38.

- Every G-L166 defining SNP tested is derived, with zero ancestral reads:
  `L166`, `L167`, `Z6516/FGC5675`, `FGC5696`, `FGC5721`, `Z6208`, `Z6219`, `Z6287`, `S19530/Z6213`.
  Five of these are transversions, so the assignment does not depend on deamination-prone sites.
- Every defining SNP of `G-Z6494` (the only child of `G-L166` in the current YFull tree) is ancestral,
  with zero derived reads: `Z6494/FGC5674` (DP 10, all MQ 60 / BQ 38), `FGC5687`, `Z6215`.
  Its subclades are ancestral too (`Z6211` DP 11, `Z6495/FGC5722` DP 10, `FT84409`).
- All ~25 markers under the provisional ISOGG `G2a2a1a2a1a1b1/b2/b2b` labels are ancestral at DP 3-13.
- Conclusion: the sample is **G-L166\*** (ISOGG-style `G2a2a1a2a1a*`).
- Evidence tables: `results/iceman_y_L166_evidence.tsv`, `results/iceman_y_deep_G_subtree.tsv`.

**Z6208 Placement Discrepancy vs Published Label (2026-07-25)**
- Nature Communications 2025 (`s41467-025-61601-8`) reports `G2a2a1a2a1a1b (G-Z6208*)`.
- That label follows the ISOGG *provisional* placement of `Z6208` (`ISOGG=G2a2a1a2a1a1b~`, note the `~`),
  which nests it inside `G2a2a1a2a1a1` = `G-Z6494`.
- Read-level data contradicts that nesting: `Z6208` is derived (DP 9) while `Z6494` is ancestral (DP 10).
  A sample cannot be inside `G-Z6494` and ancestral for its defining SNP.
- YBrowse's own YFull annotation agrees with the data: `snps_hg38.gff3` gives `Z6208` `yfull_node=G-L166`,
  i.e. `Z6208` is an L166-level SNP, not a downstream one.
- Practical rule: treat ISOGG labels carrying `~` as unranked; prefer `yfull_node` for placement.

**YBrowse ISOGG Label Reliability (2026-07-25)**
- The `ISOGG=` field in `snps_hg38.vcf.gz` is stale for parts of the G branch and produces false subclade signal.
- Observed mislabels: `FGC5696`, `FGC5721`, `FGC5687`, `Z6215` all carry `ISOGG=G2a2a1b1a` but are
  `G-L166` / `G-Z6494` defining SNPs per YFull. This is what created the spurious `G2a2a1b1a`
  derived cluster in label-based ranking.
- Do not rank on the ISOGG string alone. Resolve candidate nodes against a current tree, then test
  that node's defining SNP set explicitly.

**Can We Go Deeper Than G-L166\*? (2026-07-25)**
Going below `G-L166` requires a second individual sharing novel derived SNPs. Status of that search:

- **Sample lists on the YFull node pages could not be read reliably.** Two automated fetches of the same
  `yfull.com/tree/G-L166/` page returned conflicting sample-to-node assignments (one put
  `YF016547`/`ERS257168`/`YF096592` on `L166*`; the other put only `CGG017683` there and moved the rest
  under `G-Z6494`). The page is JS-rendered and tree indentation does not survive markdown conversion.
  **Verify node membership manually in a browser before relying on it.**
- What is consistent across reads: `G-L166*` is non-empty, and `G-Z6494` plus subclades hold ~15 samples.
- `CGG017683` (UKR) appears in this part of the tree. The `CGG` prefix is the Copenhagen GeoGenetics
  identifier scheme (Lundbeck Foundation ancient-genome panel, Allentoft et al. 2024 dense Ukraine
  sampling), so it is very likely an **ancient** sample — the highest-value comparison target found so
  far, because that data is public. Identity not yet confirmed; not resolvable via search.
- If YFull has built no subclade among the `L166*` samples, that implies they share no novel SNPs with
  each other. It does not rule out one of them sharing a branch with the Iceman, since YFull lacks him.
- Panel-based test is negative. Of 363 derived markers with DP>=4 and zero ancestral reads, every one
  resolves to the backbone `G > G-Y238 > G-P287 > G-P15 > G-PF3147 > G-PF3148 > G-PF3239 > G-L166`.
  Nothing derived below `G-L166`.
- The 82 derived markers with `yfull_node=not` are all *named backbone* SNPs (`M201`, `L31`, `L149`,
  `PF3149/3156/3237/3240`, …), not private variants. 50 of them carry `ref=Paolo Francalacci (2011)`,
  which is a nomenclature coincidence — Sardinia is G2a-rich so Francalacci named much of the G
  backbone — not evidence of shared descent.
- `ERS257168` is a weak candidate regardless: Francalacci's 1204 Sardinians are **2x low-pass**
  (only 4 at 17x). YFull's older tree (v14.04, 2018, `yfull.com/sardinians/`) placed it at `G-Z6494`;
  the current tree has it at `G-L166*`. Demotion toward the parent node is the expected signature of
  missing coverage at defining SNPs, so its `L166*` placement may be a coverage artifact rather than a
  real branch.
- The genuine high-coverage `L166*` samples are `YF016547` and `YF096592`, whose variant data is
  YFull-private. Testing them would need their novel-SNP lists.

Conclusion: `G-L166*` is the floor with currently public data. A deeper node needs a high-coverage
`G-L166*` genome to compare against.

Caveat added 2026-07-25: the "363 derived markers" test above is **label-driven** and therefore blind to
derived markers whose YBrowse ISOGG field is `unknown`/`not_listed` — 252 of them exist. See
"Blind Spot Exposed: Derived Markers With Uninformative Labels" below. The conclusion still holds for
markers ISOGG can place, but the unlabelled pool has not been resolved against a modern tree.

**CGG017683 (Goth, Ukraine) Tested Directly — Negative (2026-07-25)**
The ancient sample on the YFull `G-L166` page resolves to public data, and was tested read-level.

- Identity: `CGG017683`, Ukraine, ~550 CE (300-800 CE), Migration Period, Gothic cultural context,
  mtDNA `V1a10`.
- Data: ENA sample `SAMEA117822726`, run `ERR14752008`, study `PRJEB87274`
  ("Tracing the Spread of Germanic Languages using Ancient Genomics", Centre for GeoGenetics,
  712 ancient genomes). **First public 2025-04-08** — earlier searches predate release.
- Submitted BAM is GRCh37/hs37d5-aligned (`1..22,X,Y,MT`) and ships with a `.bai`, so
  `samtools view -b <https-url> Y` streams the chrY slice without downloading all 27.3M reads.
- **Coverage is the limiting factor**: 54,685 chrY reads, 2,413,071 bases = **0.041x on chrY**.
  P(any read at a given site) = 4%. All reads pre-filtered to MQ>=25.
- Marker test (56 deep-G markers lifted hg38->hg19): 11 covered, all single-read. Derived at
  `S19530/Z6213` (an L166-level SNP, consistent with the YFull assignment) and at `Z6504`
  (single deamination-prone read where the Iceman is ancestral at DP 7 — treat as damage).
- Shared-variant test: of the Iceman's 438 high-quality chrY derived sites (PASS, `GT=1/1`,
  DP>=5, GQ>=30), 55 are covered in `CGG017683` and 51 share the allele — but **all 51 are
  catalogued YBrowse markers**, i.e. upstream/L166-level SNPs both carry by common ancestry.
  **Zero novel shared variants.**
- The Iceman has only ~20 novel (uncatalogued) HQ chrY sites (17 by the VCF panel, 20 by the fuller
  GFF3 index). `CGG017683` covers **3** of them and is **ancestral at all three**
  (`hg38 7899558`, `hg38 11414525`, `hg38 19647870`).
- Verdict: no shared sub-L166 branch detected. Direction of evidence is right but power is very low
  (3 informative sites, 1 read each). This does not refute a shared branch, it fails to find one.

  *(Updated 2026-07-25: this bullet pair originally said 2 sites / "ancestral at both", from the
  pre-MAPQ candidate list. The MAPQ-filtered re-derivation later in this file re-tested `CGG017683`
  against the surviving candidates and found 3 covered, ancestral at all 3. The verdict is unchanged
  and the power is still very low; only the count moved. The later section is authoritative.)*

**chrY DoC Denominator Gotcha (2026-07-25)**
`CGG017683` reports **0.192x** chrY in the study but **0.041x** if computed over the full 59.37 Mb contig.
Both are the same 2,413,071 mapped bases; the study divides by a ~12.57 Mb callable/MSY denominator
(2413071 / 12.57e6 = 0.192). Y heterochromatin has no reads, so the full-contig figure understates real
site coverage by ~4.7x. **Always state the denominator.** Empirical hit rates confirm the callable
figure: 19.6% of G-branch markers and 12.8% of the Iceman's variant sites were covered, consistent with
`1 - e^-0.192 = 17.5%`, not 4%.

**PRJEB87274 Sibling Samples Screened — All Dead Ends (2026-07-25)**
Same Crimean site, same study as `CGG017683`. Screened via remote indexed chrY fetch + hg38->hg19
liftover + pileup (~2 min per sample, no bulk download):

| sample | total bases | chrY DoC (12.57 Mb) | verdict |
|---|---|---|---|
| `CGG017681` | 6.12 Gbp | 0.008x | **female** — 2,251 chrY reads despite 5x the data |
| `CGG017682` | 3.51 Gbp | **0.623x** | male, good depth, but **not this lineage** |
| `CGG017683` | 1.18 Gbp | 0.192x | only G-L166 one, too shallow to test |

`CGG017682` is ancestral at `L166`, `L167`, `L91`, `PF3148`, `Z6043` — reliable transversion/clean calls.
Its apparent derived G-panel hits are **all single-read `C>T`/`G>A`**, i.e. deamination artifacts; `U1` is
ancestral at 2 reads. Useful validation that the damage down-weighting in `y_path_rank.py` targets a real
failure mode.

**Method Positive Control (2026-07-25)**
The paper assigns `CGG017682` to `M173 > L146 > M459 > M417 > PF6162 > S202 > Z94 > Z2124 > S4576 >
Y57 > YP1269` (R1a-Z2124). Re-derived independently here from the public BAM via remote indexed fetch +
hg38->hg19 liftover + pileup: 5 of 11 path markers covered, **4 of 4 informative ones DERIVED** —
`L146/M420` (T>A transversion, 2 reads), `M459`, `Z94`, and terminal `YP1269`. `M173` uncovered.
This validates the whole remote-screen chain against a published answer, and confirms that the sample's
apparent "derived G" hits were deamination artifacts: it is R1a, not G.

**PRJEB87274 Search Closed (2026-07-25)**
Per the study supplement, `CGG017683` is the **only G2a2a individual among all 712 genomes**. There is no
second candidate in this dataset, so a systematic BAM screen is unnecessary — the one G2a2a sample has
already been tested and is depth-limited at 0.192x. Search within this study is exhausted.

**Independent Support for the L166-Terminal Correction**
The study assigns `CGG017683` the path
`P15 > M3308 > PF3147 > PF3148 > PF3177 > L91 > Z6488 > PF3239 > L166`.
Two things follow:
- The Iceman is **derived at all nine**, with **zero ancestral reads at any of them**:
  `P15` DP16, `M3308` DP10, `PF3147` DP2, `PF3148` DP8, `PF3177` DP11, `L91` DP7,
  `Z6488` (chrY:22225244) DP3, `PF3239` DP10, `L166` DP3. The two individuals traverse an identical
  path and terminate on the same node.
- The Copenhagen pipeline **terminates that sample at `L166`** — it does not invoke `Z6208` and does not
  place the sample inside `Z6494`. A current aDNA pipeline handling a G-L166 genome stops exactly where
  our read-level analysis says the Iceman stops, which is independent support for reading the published
  `G2a2a1a2a1a1b (G-Z6208*)` label as an inherited ISOGG provisional artifact.

**Where a Branch-Mate Might Actually Be Found**
`G-L166` is an Anatolian-farmer-derived lineage, so Neolithic Europe is a far better hunting ground than
a Migration Period Germanic dataset. Candidate pools:
- AADR (Allen Ancient DNA Resource) — curated compendium with Y calls.
- aYChr-DB — ~1,797 ancient Y-chromosome samples.
- Published Neolithic G2a series (LBK, Anatolian farmers, Iberian/French/Balkan Neolithic), which are
  G2a-rich by comparison.
Requirement: enough chrY depth to cover the 8 usable candidates in
`results/iceman_y_novel_branch_candidates.tsv` — the MAPQ-filtered set described under
"Novel-Branch Candidates Re-derived With MAPQ Filtering", of which only 1 is a transversion. The
site list is committed as `markers/iceman_novel_candidates_usable8.tsv`. At `CGG017683`-level 0.192x the expected
yield is ~1-2 informative sites, so target samples above ~1x chrY.

**FTDNA Block Tree Cross-Check (2026-07-25)**
FTDNA's Block Tree groups phylogenetically-equivalent SNPs into blocks, and its blocks around `Z6488`
are *finer* than ISOGG's, which lumps all of them at `G2a2a1a2a`. Two blocks were tested at read level:

| FTDNA block | marker | chrY hg38 | anc>der | DP | anc | der | call |
|---|---|---|---|---|---|---|---|
| `G-Z6488` | Z6488/FGC7739 | 22225244 | T>C | 3 | 0 | 3 | DERIVED |
| `G-Z6488` | FGC2271/Z6484/S25082 | 20907600 | G>A | 6 | 0 | 6 | DERIVED |
| `G-Z6488` | PF3237 | 14905951 | G>A | 4 | 0 | 4 | DERIVED |
| `G-PF3238` | PF3238 | 15144551 | G>A | 5 | 0 | 5 | DERIVED |
| `G-PF3238` | FGC2274 | 5094513 | G>A | 1 | 0 | 1 | `nocall_damage_prone_1read` (was "DERIVED (weak)") |
| `G-PF3238` | Y232445 | 11361597 | C>T | 9 | 0 | 9 | DERIVED |
| `G-PF3238` | Z6118 | 3824595 | T>C | 2 | 0 | 2 | DERIVED |
| `G-PF3238` | Z6119 | 3824608 | A>G | 3 | 0 | 3 | DERIVED |

8/8 derived, 0 ancestral reads anywhere. The Iceman passes cleanly through both FTDNA blocks. This does
**not** by itself go deeper than `G-L166` — `PF3238`/`PF3239` are consecutive PF numbers and the block
most likely sits *between* `G-Z6488` and `G-L166` (ISOGG `G2a2a1a2a1`). Confirming that requires reading
the parent/child order off FTDNA's Block Tree directly; the site is not machine-accessible from here.

Note `Z6118`/`Z6119` are 13 bp apart and can share a read — treat them as one observation, not two.

**Blind Spot Exposed: Derived Markers With Uninformative Labels (2026-07-25)**
The `G-PF3238` block markers carry YBrowse `ISOGG="unknown"` / `"G_(not_listed)"`, so **both** prior
searches missed them:
- the G-subtree pileup selected markers by ISOGG `G*` label — these have no usable label;
- the private-branch scan kept only variants absent from the panel — these are *in* the panel.

Systematic re-count over the 835 derived panel markers passing GQ>=20/DP>=3: **125 carry an
uninformative ISOGG label on *every* alias**, 107 of them at sane depth (DP 4-20). Written to
`results/iceman_y_unlabelled_derived_markers.tsv` with hg38 + hg19 coordinates, DP, GQ, and the raw
YBrowse `HG`/`ISOGG` strings.

Count definition matters: a first pass flagged 252, but that counted any marker with *at least one*
unplaceable alias. A marker like `Z6284,S17295` labelled `G_(not_listed),not_listed` is genuinely
unplaceable; one labelled `G2a2a1a2a1,unknown` is not. The strict rule (all aliases unplaceable) gives
125 and is what `is_unplaceable_label()` in `y_haplo_from_markers.py` implements.

This is the pool in which a below-`L166` marker would hide, and FTDNA's Block Tree / YFull are exactly
the resources that place such markers. Depth filter matters: the top of the list by DP
(`A2472` DP 906, `BY49567` DP 264, `BY31542` DP 113, `A2500` DP 103) is the same collapsed-repeat
artifact class documented under "Marker-Level Artifacts in Repetitive chrY" — reject anything far above
the ~8x chrY mean. Positional clusters at ~11.5 Mb and ~26.6 Mb are likewise suspect.

**Tooling Added to Close the Blind Spot (2026-07-25)**
- `annotate/build_marker_index.sh` — collapses the 796 MB / 3.1M-line YBrowse GFF3 into
  `resources/marker_index.tsv.gz` (name, chrom, pos, anc, der, isogg, ycc, yfull_node, ref, comment).
  Takes ~2 min to build. This is the actual root cause of the label-driven habit: without a name index,
  every by-name lookup was a full scan of an 800 MB file, so selecting by ISOGG label was the only
  practical option.
- `annotate/y_markers_pileup.py` — takes marker *names* (`--markers` or `--marker-file`), resolves them
  through the index, and reports ancestral/derived read counts straight from the BAM with **no label
  filtering**. Reports `not_in_catalogue` explicitly rather than silently dropping. A block of SNP names
  pasted from FTDNA's Block Tree can now be tested in seconds:
  ```
  annotate/y_markers_pileup.py --bam "$BAM" --ref "$REF" \
      --marker-file block.txt --label G-PF3239 --out results/block.tsv
  ```
  It extracts reads by BED first (see the ops gotcha in RUNLOG) so a 24-marker block runs in ~6 s.
- `y_haplo_from_markers.py` now emits a fourth output `<prefix>.unplaceable_derived.tsv` and a
  `derived_unplaceable` summary count, so this pool can never again be invisible in a normal run.

Regression check: re-running every hand-built pileup from this session through
`y_markers_pileup.py` reproduced all depths and calls exactly. Consolidated in
`results/iceman_y_ftdna_block_evidence.tsv`.

**FTDNA Block Walk: G-PF3239 vs G-FGC2315 (2026-07-25)**
All 24 names in FTDNA's `G-PF3239` block resolved in the catalogue — **22 DERIVED, 0 ancestral**, plus
two `mixed` sites that are clearly paralogous (`Z6309` chrY:11265625, 9 anc / 8 der at DP 17 on a
haploid chromosome; `Z6312` chrY:11436604, 2/4). Five of the derived calls are transversions
(`FGC5664`, `FGC5676`, `FGC5681`, `Z6126`, `Z6492`), so the block is damage-independent.

The sibling block is excluded: `FGC2315` (chrY:7679120, T>C) is **6/6 ancestral**.

Independent sanity check from ISOGG's own numbering: `PF3239` is `G2a2a1a2a1` and `FGC2315` is
`G2a2a1a2a2` — siblings. The Iceman is derived on the `...a1` side and ancestral on the `...a2` side,
exactly as the FTDNA split predicts.

Six markers in this block carry ISOGG `G2a2a1b1` (`FGC5676`, `FGC5681`, `Z6126`, `Z6129`, `Z6283`,
`Z6486`) while their `yfull_node` is `G-PF3239` — another instance of ISOGG label staleness. **Prefer
`yfull_node` over `isogg_haplogroup` when the two disagree.**

**Still upstream of L166.** `PF3239` = ISOGG `G2a2a1a2a1`; `L166` = `G2a2a1a2a1a`, i.e. the next level
down. So the walk `G-Z6488 > G-PF3238 > G-Z6128 > G-PF3239` is still descending *toward* `L166`, not
past it. The split that matters is whatever `G-L166` itself divides into.

**Catalogue-Wide YFull Node Census (2026-07-25) — supersedes the label-driven test**
With `resources/marker_index.tsv.gz` in place it is cheap to ask the inverse question directly: over
*every* marker the Iceman is derived at (834 HQ positions), which YFull nodes do they belong to? This
replaces the earlier ISOGG-label test, which could only see markers ISOGG had ranked.

Every G-node represented, by marker count:
`G-PF3147` 34, `G-Y238` 31, `G-P15` 26, `G-PF3239` 23, `G-P287` 17, `G-PF3148` 12, **`G-L166` 12**,
`G-Z6484` 4, `G-PF3233` 3, `G-L1259` 2, `G-Z6284` 1, `G-Z6128` 1, `G-L91` 1, `G-PF3177` 1.

All are on the backbone at or above `G-L166`. **Nothing below it.** Two apparent exceptions were chased
down and both are join artifacts, not deeper nodes:
- `G-Z12221` (`FGC2306`, chrY:6777446, A>T) shares a position with `Z6128` (A>C). The Iceman carries C,
  so he is `Z6128`-derived and *not* `FGC2306`-derived. Multiallelic site; a position-keyed join
  attributes it wrongly.
- `G-Y171774` (chrY:15732901, C>A) — chrY:15732901 carries **four** markers with differing alleles.
  The Iceman's C>G matches `MF792575`, whose ISOGG is `O` and whose `yfull_node` is `not`, i.e. a
  recurrent mutation in an unrelated lineage. Not informative for G.

Lesson: when joining variants to the catalogue **by position, always match the allele too**, and expect
multiple markers per position with conflicting anc/der.

**Third-Allele Sites (a third blind spot, checked and closed)**
Sites where the Iceman carries an allele that is neither the catalogued ancestral nor any catalogued
derived allele are invisible to both prior scans — the novel-variant scan skips them (position *is*
catalogued) and the marker scan skips them (allele does not match). Catalogue-wide there are **6**.

First-pass filtering found 50, but 44 were the familiar **hg38-chrY-reference-carries-derived-alleles**
artifact: the Iceman shows as "variant" at `M269`, `P312`, `U152`, `M173`, `M9` etc. purely by being
ancestral where the reference is derived. Filter against the catalogued **ancestral** allele as well as
the derived one, or this swamps the result.

Of the 6 real ones: 3 sit at chrY:56.70-56.75 Mb and 1 at 26.66 Mb (known collapsed-repeat zones,
reject); chrY:10829384 A>T is rejected on mapping quality (4 reads at MAPQ 0, nothing above MAPQ 40 --
same failure mode as `Z6519`); chrY:20327833 G>T survives at DP 5 with 4 reads MAPQ 60, a transversion,
but **all 5 reads are forward-strand**, so treat as marginal. Added to the novel-branch candidate set
with that caveat.

**Conclusion:** `G-L166*` is the floor *within the catalogue*. Every catalogued marker, placeable or
not, has been tested. What the census cannot speak to is markers the catalogue does not contain — see
next section.

**Catalogue Refresh + Exhaustive G-L166 Subtree Test (2026-07-25) — the definitive negative**
Refreshed YBrowse from upstream (`snps_hg38.gff3`, 2026-07-24, 817,252,415 bytes vs the local
2026-02-23 copy at 796,023,961). The rebuilt index has **3,206,354 rows, +84,355 on February**.

Three results, all negative and all worth having:
1. **None of the 11 surviving novel-branch candidate positions has acquired an SNP name.** The
   uncatalogued set is byte-identical against the fresh catalogue: still the same 20 positions.
2. **The `G-L166` subtree has not changed at all** in five months — same nodes, same marker counts.
3. **The Iceman is ancestral across the entire subtree.** All 93 markers under `G-L166` were tested
   (nodes `G-Z6494`, `G-Z6036`, `G-Z6211`, `G-Z6212`, `G-Z6214`, `G-Z6285`, `G-FGC5697`,
   `G-CTS12768`): **88 ancestral, 2 no-coverage, 0 credible derived.**

The three non-ancestral calls are all artifacts:
- `FT172194` / `Y126330` are the *same position* (chrY:2891795), DP 1, C>T deamination-prone, and
  3 of 4 reads there are MAPQ 0. A single damage-pattern read on bad mapping.
- `Y38112` (chrY:13029000) is `mixed` at 2 anc / 4 der with 3 of 11 reads MAPQ 0 — paralogous.

Many of the ancestral calls are transversions at DP 7-16, so this is not a coverage-limited result.

**This is what the FTDNA Block Tree walk was reaching for, obtained from the catalogue instead.**
`G-L166` does have children — eight of them, all well populated with modern testers — and the Iceman
belongs to none. `G-L166*` is therefore a genuine terminal placement, not an artifact of missing data
or missing annotation.

**FTDNA G-PF3147 Project Page Tested (2026-07-25)**
A saved FTDNA public *group project* page for `G-PF3147` / `G-L166` was tested against the Iceman.
25 member rows, 3,788 distinct SNP names, 18 distinct terminal haplogroups (all `G-BY*`).

Result: **all 18 terminal SNPs are ancestral in the Iceman** — several at good depth and on
transversions (`BY220097` DP 17 A>T, `BY202787` DP 16, `BY66658` DP 13, `BY106602` DP 12 G>T,
`BY111300` DP 12 G>C). Across all 3,788 names there is **not one derived call carrying a G-lineage
node off his known backbone**. The 1,733 derived hits are backbone or `A0-T`-level markers that every
non-A00 human carries.

**The page is filtered to G-L166-positive members only** (`G-PF3147` is the project; `G-L166` is the
filter). So this is not a scattered G2a sample — it is 25 confirmed `L166+` men spread across 18
terminal branches, and the Iceman is ancestral at all 18. Corroborating detail: members show `Z6211+`
and `FGC5697+`, both `L166` children from the subtree test, and their terminal nodes (`G-Y97527`,
`G-Y140820`, `G-BY123198`, `G-Y53493`, `G-Y154380`, `G-Y256820`, …) are *finer* than anything
YBrowse's `yfull_node` field records. Do not read the `Confirmed SNPs` column as the filter — it
spells out `L166` for only 5 of 25 rows and is selective display.

This is the best-powered living-carrier test available: every sampled branch of `G-L166` that has
modern representatives, and the Iceman is on none of them. Caveat: public project entries only, so it
is a subset of FTDNA's L166 testers.

Saved: `results/iceman_y_ftdna_project_terminals.tsv` (SNP names + the Iceman's own calls only).

**Handling note — the source page must not be committed.** These pages carry `Kit Number`,
`Last Name`, and `Paternal Ancestor Name` columns for living project members. The file was moved to
`<private_refs_dir>` (outside the repo) and `.gitignore` now guards `*.html`, `*.htm`,
`FamilyTreeDNA*`, `YFull*`. Only aggregate SNP names and the Iceman's own genotypes are in-repo, which
is consistent with the Data Policy above. Note the repo had **no** prior guard for saved web pages —
the file was untracked but committable.

**Catalogue Completeness: the VCF and GFF3 Disagree (2026-07-25)**
`resources/snps_hg38.vcf.gz` (2,920,250 records) and `resources/snps_hg38.gff3` (3,121,999 rows,
2,738,775 distinct chrY positions) are **not the same marker set**. The earlier novel-variant scan ran
against the VCF and called 17 positions novel; re-running against the GFF3-derived index gives a
different set of **20**. Ten of the original 17 (the chrY:56.69-56.88 Mb block) are in fact catalogued
in the GFF3, while 13 positions absent from the VCF had never been evaluated at all.

Treat the GFF3 as authoritative; the VCF has been through `sanitize_marker_vcf.sh` and is lossy.
Anything claiming a variant is "novel" must say which catalogue it was checked against.

Both files are from the 2026-02-23 refresh — roughly five months stale as of this work. FTDNA registers
new SNP names (`FT*`, `BY*`, `MF*`) with YBrowse continuously, so `annotate/fetch_ybrowse_markers.sh`
should be re-run before any future claim that a variant is uncatalogued.

**DeepVariant Writes No Provenance, So The Filename Carried All Of It (2026-07-26)**
The complete non-boilerplate header of both A/B outputs is:

```
##fileformat=VCFv4.2
##DeepVariant_version=1.10.0
```

Version and nothing else — no command line, no `--ref`, no `--reads`, no `--disable_small_model`.
**The two branches are indistinguishable at the header level.** The filename was carrying the
entire provenance, and the filenames did not survive:

| RUNLOG says | on disk | index | written |
|---|---|---|---|
| `iceman_sm_on.vcf.gz` | `~/iceman.vcf` (BGZF content, `.vcf` name) | none | Feb 18 |
| `iceman_sm_off.vcf.gz` | `~/iceman-nosmall.vcf.gz` | `.tbi` | Feb 21 |

Neither recorded name survived, and the Feb 18 pair carries compressed content under a bare
`.vcf`/`.gvcf` extension, so `cat`/`grep` return binary. Identification was only possible
**forensically, from the data** — `PASS,GT=1/1,DP>=5,GQ>=30` returning exactly 438 on one branch
and 0 on the other, and chrY `GQ` max 58 vs 25. That is a reconstruction, not a provenance trail,
and it worked only because a distinguishing statistic happened to exist. Had the two branches had
similar GQ ranges these files would now be permanently unidentifiable.

**Design conclusion: content sniffing is necessary but not sufficient; use hash-as-identity.**
A tool that does not record its own parameters cannot be made self-describing after the fact, so
the wrapper must do it: emit a sidecar manifest per run carrying the exact command, container
digest, input checksums, output `sha256` and timestamp. An inventory tool then matches files to
manifests **by checksum, not by name**, after which renaming is harmless — the hash is the
identity and the name is only a label.

For legacy files with no manifest, record the *distinguishing statistic* alongside the
identification so it is re-derivable rather than a one-off deduction. For these two that is
chrY hom-alt `GQ` max: **58 = small model on, 25 = off**.

**`GQ>=30` On chrY Is A Small-Model Artefact (2026-07-26)**
The DeepVariant A/B outputs were located at `~/iceman.vcf` (small model on) and
`~/iceman-nosmall.vcf.gz` (off), both DeepVariant 1.10.0, 33,546 chrY records each. Note
`iceman.vcf` and `iceman.gvcf` are **gzip/BGZF despite the extension** — `cat` yields binary.

Hom-alt chrY SNVs, all `PASS` in both branches:

| branch | GQ median | GQ **max** | sites at `GQ>=30` |
|---|---|---|---|
| small model on | 17 | **58** | **438** |
| small model off | 17 | **25** | **0** |

**Without the small model, GQ never exceeds 25 on this data.** `GQ>=30` is therefore not a neutral
quality bar; it is unreachable unless the small-model path ran. The medians are identical, so the
difference is entirely in the confident tail — consistent with a fast heuristic path emitting
high confidence on "easy" sites, which is exactly the behaviour whose validity on damaged aDNA is
in question.

This identifies the provenance of the existing analysis: the filter `PASS, GT=1/1, DP>=5, GQ>=30`
returns **exactly 438** on the small-model branch, matching the count recorded above, and **0** on
the other. So `~/iceman.vcf` is the file the earlier chrY work used.

**Scope of the caveat.** It touches the novel-variant discovery and the 438-site `CGG017683`
comparison. It does **not** touch the `G-L166` result, which comes from `y_markers_pileup.py`
reading the BAM directly — read counts and MAPQ, no genotyper and no GQ anywhere in the chain.
The headline finding is VCF-independent and unaffected.

Only **7 of the 21** novel candidates appear in the 438 set, because several sit at DP 4 and the
novel scan evidently used a looser filter than the shared-variant analysis. The exact filter that
produced the 21 is still not recovered — the reason `annotate/y_novel_scan.py` should derive
candidates from the BAM rather than from a VCF, where the quality scale is model-dependent.

**THE COMPILER HEDGED; THE PRESENTATION LAYER DROPPED IT (2026-07-26)**
The single most important thing in this whole thread. In `all-ancient-dna.txt`, of the 22 samples
under `Culture_Grouping = "Horgen culture"` called `G-L166(*)` in `Y-Haplotree-Variant`, **21 carry
a parenthesised ISOGG node**:

    ISOGG2019 = G2a2a1a2a1(a)      # G2a2a1a2a1 = PF3239, G2a2a1a2a1a = L166

The parenthesis means *PF3239, possibly L166* — provisional. The **one** unhedged entry is
`TU876` (`SX10`) at `G2a2a1a2a1a`, and that is precisely the sample Furtwängler 2020 genuinely
reports as `L166`. The compiler distinguished the confirmed case from the 21 provisional ones,
sample by sample, and got it right.

`Y-Haplotree-Variant` then flattens both forms to `G-L166*` / `G-L166`. That is the column the web
front-ends expose (`haplotree.info` searches on it) and the one downstream readers quote.

**So the failure is not fabrication and not the compiler's.** The uncertainty was recorded and is
still in the file; it is lost one column over, at the presentation layer. Every claim built on
`Y-Haplotree-Variant` inherits 21 provisional calls as though they were 21 confirmed ones.

This is the same structural error as the Iceman's `G-Z6208` label, arrived at from the opposite
direction: there, a real terminal-derived SNP was read as a terminal placement; here, a hedged
node is read as a confirmed one. Both times the evidence survives and the qualifier does not.

Practical rule for this repository: **never cite `Y-Haplotree-Variant` without checking
`ISOGG2019` for parentheses.**

**Why the hedge is dropped: a documented naming policy selects the column that cannot carry it.**
`github.com/Dimetrodon2026/Archaeogenetics` (README only, plus a "Links to G25 collections"
folder; no code, no pipeline) states the convention used for the G25 collections: **YFull naming
for haplogroup G**, FTDNA naming for most other lineages, theytree for some East Asian branches.

In `all-ancient-dna.txt`, across the 29 Horgen-grouped males with a Y call:

    parenthesised hedge present in ISOGG2019 : 21
    parenthesised hedge present in YFull     :  0

The hedge is not merely absent from `YFull` — it is *inexpressible* there. ISOGG longhand can
write `G2a2a1a2a1(a)` = "PF3239, possibly L166". YFull shorthand asserts the node and qualifies
only what lies below it. Converting ISOGG → YFull is lossy in exactly the dimension under dispute,
and the policy is a perfectly defensible one that happens to land on the lossy side.

**The asterisk then inverts the meaning.** Standard YFull/ISOGG convention: `*` marks a sample in
a clade but in none of its known subclades — a *terminal* placement, which is precisely the
Iceman's real status (`G-L166*`). In this file the relationship is reversed:

| | ISOGG2019 | YFull |
|---|---|---|
| `SX10` — publication reports `L166` | `G2a2a1a2a1a` | `G-L166` |
| the other 21 — provisional | `G2a2a1a2a1(a)` | `G-L166*` |

The confirmed sample carries the bare name; the provisional ones carry the asterisk. A reader
applying the standard convention therefore reads the file backwards, seeing `G-L166*` and
understanding "confirmed terminal L166, like Ötzi" — which is the circulating claim verbatim.

*This asterisk reading is an **inference** from the `SX10` contrast. No published legend for the
file has been located. Do not present it as documented.*

**On the "Dimetrodon → FTDNA/GenArchivist → AADR" provenance line** attached to the circulating
graphs ("deduped per individual, dates verified in-window. 344 males / 11 groups"), reported by
the user 2026-07-26; the post itself is unreadable here (x.com returns HTTP 402 unauthenticated):

- The QC claims are **sample bookkeeping**, not haplogroup validation. Deduplicating individuals
  and checking radiocarbon dates fall inside a culture's window are both real and worth doing, and
  neither says anything about how `G-L166` was called for any skeleton.
- `AADR` plausibly supplies samples, dates and individual identity — the auditable part.
  `FTDNA/GenArchivist` is the layer supplying Y calls, and GenArchivist is exactly where the
  "analyses of BAM files by hobbyists" named in `haplogroup.info`'s own methodology are posted.
- **Not verified:** what Y-haplogroup resolution AADR's `.anno` actually carries. The AADR
  Scientific Data paper does not document its Y columns, and the file has not been downloaded.
  Do not assert that AADR cannot be the label source until that is checked.

**SOURCE OF THE "7 OF 10 ARE G-L166" CLAIM IDENTIFIED (2026-07-25)**
It is the **"All Ancient DNA" compilation** at `haplogroup.info` — distributed via
`indo-european.eu/ancient-dna/`, and the backend of `haplotree.info` (same `Y-Haplotree-Variant`
column name). 16,972 rows; retrieved `haplogroup.info/all-ancient-dna.txt` 2026-07-25.

Unlike its web front-ends it *does* cite per-sample sources. The chain:

1. `Source = FurtwaenglerNatCommun2020` — same paper we have been reading.
2. The Y calls are **not the paper's**. The compilation re-derives them; the page states
   haplogroup inference by **Yleaf v2** plus "analyses of BAM files by hobbyists and online
   informal reports of research papers in preparation". Each row credits a
   `Responsible-SNP` — for the whole Swiss series, **Milan Rajevac**.
3. That re-derivation relabels **6 of the 7** published-`PF3239` Oberbipp males as **`G-L166*`**
   (`MX187`, `MX209`, `MX210`, `MX211`, `MX212`, `MX213`). Only `MX182` stays `PF3239` — and its
   `NRY` count is **2**, against 38,316 for `MX210`. `MX183`, published `Z6488`, is moved *deeper*
   to `PF3239`.
4. It assigns `Culture_Grouping = "Horgen culture"` to **Aesch, Oberbipp, Muttenz and Rapperswil
   alike** — three of which the publication does not label Horgen at all.

Querying that grouping reproduces the claim's shape exactly: **49 samples, 29 with a Y call, of
which 22 are `G-L166(*)`** — including 12 of 14 at Aesch. Restricted to Oberbipp it is 6 of 10 in
the current file; the 2022 snapshot on `indo-european.eu` is the likely origin of "7 of 10", and
`SX10` (`G-L166`, Rapperswil) is a candidate seventh under this grouping.

So the count is real, the culture grouping is the compilation's, and the `L166` label is a
hobbyist re-call. The publication says `PF3239`. Nothing was fabricated; a re-derived label was
read as if it were the published one — the same conflation this project documents for the Iceman,
one layer further out.

**The compilation's own kinship note undercuts the population reading.** Verbatim from `MX211`:

> Switzerland, Oberbipp Horgen **Family A (9 members)** (MX150-MX187 have a father-son
> relationship (order unknown). MX187-MX212 have a father-son relationship (order unknown). …
> MX183-MX211 are brothers … MX209 (who is a 1st degree relative of MX219)

Five of the seven relabelled men (`MX182`, `MX187`, `MX209`, `MX211`, `MX212`) sit inside one
named nine-member family; only `MX210` and `MX213` fall outside it. Even taking the `L166` calls at
face value, this is largely one patriline, not a population.

**And our read-level test disagrees with it on the one sample tested so far.** The compilation
calls `MX210` `G-L166*`. Our mapping of the raw FASTQ finds it **ancestral at `L166` (3/3 reads)
and `L167` (2/2)**, both transversions at 0% MQ0 — while derived at `Z6219` (4/4). That is
recorded below as `CONFLICT`, not as a refutation, and it is one sample of nine. But it is the
first direct read-level check of a `G-L166*` call from this compilation, and it did not confirm it.

**`L166` And `L167` Are Recurrent Sites, And MX210 Conflicts (2026-07-25) — interim, 5 of 9 mapped**
The marker index carries multi-haplogroup YFull assignments for the two SNPs the node is named for:

| marker | chrY hg38 | class | isogg | yfull_node |
|---|---|---|---|---|
| `L166` | 21843737 | transversion | G2a2a1a2a1a | **`G-L166&J-Y29712`** |
| `L167` | 21843756 | transversion | G2a2a1a2a1a | **`G-L166&I-Y92994&A00`** |
| `Z6219` | 13782251 | transition | G2a2a1a2a1a | `G-L166` |

`L166` and `L167` are therefore **recurrent** — the same positions define nodes in haplogroups J,
I and A00. That is a catalogue-level caveat on the two markers carrying the node's name, and it
was not previously recorded here. It does not weaken the Iceman result, which is derived at all
nine including six non-recurrent ones, but it matters wherever those two carry a call alone.

Note also the set is not internally consistent in ISOGG terms: `FGC5696` and `FGC5721` are
`G2a2a1b1a`, a different path from `G2a2a1a2a1a`; `Z6208` is `G2a2a1a2a1a1b~`, one node deeper;
`S19530` is `not listed`. All nine share `yfull=G-L166`. YFull groups where ISOGG scatters.

**`MX210` (0.1434x, best of the set) returns CONFLICT**, and it is not obviously artefactual:
- `Z6219` C>T: DP 4, **4 derived / 0 ancestral**, 0% MQ0, 4/4 at ceiling.
- `L166` C>A transversion: DP 3, **3 ancestral / 0 derived**, 0% MQ0, 3/3 at ceiling.
- `L167` T>A transversion: DP 2, **2 ancestral / 0 derived** — but its 2 reads are a *subset* of
  L166's 3 (reads start at 21843678/68bp, 21843708/56bp, 21843721/48bp; only the latter two reach
  21843756). So this is 3 independent molecules, not 5.
- Six of nine sites uncovered. `PF3239`, the positive control, uncovered; `P15`, `M3308`, `L91`,
  `Z6043` are derived, so the backbone is partly confirmed.

Deamination was tested for and does not explain `Z6219`: of its 4 reads, 3 carry the T well
inside the molecule (offsets 30, 28, 15 from the 5' end; 37, 15, 34 from the 3'), and only one
sits within 5 bp of a terminus. Under UDG-half that is not a damage pattern.

So the two damage-immune transversions say ancestral on 3 molecules, and one non-recurrent
transition says derived on 4. **No call is being made from this.** `CONFLICT` is the registered
verdict for a sample with derived and ancestral calls in the same set, and it stands until the
remaining four samples are mapped and the other six sites have coverage somewhere.

**The MAPQ-60 Test Is Meaningless For `bwa aln` Data (2026-07-25)**
`bwa aln` caps MAPQ at **37**; `bwa mem` uses **60**. Measured, not assumed: of 30,825 chrY reads
in `MX210`, 30,625 sit at exactly 37 and the file maximum is 37, while the Iceman `bwa mem` BAM
reaches 60 in the same window.

The `n_mq60 == 0 -> MARGINAL` rule in `site_qc()` was calibrated on the Iceman BAM and therefore
evaluated to "no uniquely-mapping reads" for **every read of every site** in the Swiss capture
data. It flagged 100% of sites, including ones at 0% MQ0 with full-ceiling mapping — a quality
test no `bwa aln` alignment can ever pass. Nothing errored; sites silently failed a bar that did
not apply to them.

Fixed by detecting the aligner from `@PG` (`bwa samse`/`aln` -> 37, `bwa mem` -> 60) and counting
reads at *that* ceiling. Columns `n_mq60` -> `mq_top` + `n_mq_top`, so the threshold in force is
recorded per row rather than implied. Iceman re-validation is unchanged: ceiling detected as 60,
**21/21 verdicts, 0 field mismatches**.

Two related gaps closed at the same time:
- `y_markers_pileup.py` emitted **no MAPQ columns at all**, although the prereg calls the audit
  mandatory and it is the tool that runs the primary L166 test. It now emits the same audit.
  (Strand counts stay unavailable there — its `parse_bases()` is strand-agnostic — so the
  single-strand test is skipped rather than reported wrongly.)
- Call rules, mutation classing, MAPQ handling and region flags now live in `annotate/ylib.py`
  and are imported by both tools, because they had already drifted apart once.

**Ancient `G-L166` Candidate Pool From haplotree.info (2026-07-25)**
`https://haplotree.info/maps/ancient_dna/slideshow_samples.php?searchcolumn=Y_Haplotree_Variant&searchfor=G-L166&ybp=500000,0`
(dataset "All Ancient DNA v. 2.07.26"). **7 samples**, 0.05% of their corpus:

| sample | ybp | site | culture |
|---|---|---|---|
| `I5118` | 5100 | Mezőcsát-Hörcsögös, Hungary | Hungary_Baden LCA, Viss group |
| `TU876` (= `SX10`) | 4550 | Rapperswil Zürichstrasse, CH | Horgen *(their attribution)* |
| `I14677` | 4329 | Serra Crabiles, Sardinia | Sardinia_Bell_Beaker |
| `UNTA58_68Sk1` | 4302 | Haunstetten, Unterer Talweg 58–62, DE | Bell Beaker Lech Valley |
| `I15942` | 4297 | Anghelu Ruju, Sardinia | Sardinia_Bell_Beaker |
| `E09538` | 4287 | Unterer Talweg 58–62, Augsburg, DE | BK_Germany_BAV |
| `I14678` | 4247 | Serra Crabiles, Sardinia | Sardinia_Bell_Beaker |

This answers the "candidate pools" open question above with named samples. `I5118` at 5100 ybp is
near-contemporary with the Iceman. `UNTA58_68Sk1` and `E09538` are the same site at nearly the same
date and may be one individual under two aggregation IDs — check before counting them as two.

**`CGG017683` is absent from this list.** It was reachable through YFull, which is where the
impression that it was the only close ancient `L166` came from. YFull is *selective* — it ingests
submitted BAMs meeting its own criteria — whereas haplotree.info draws on AADR-style compilations.
Neither is a superset of the other, so "the only one on YFull" was never evidence of "the only
one". Any future claim about how many ancient `L166` exist has to name the corpus it counted.

*Correction (2026-07-26): the 7 above are the **unhedged** `G-L166` rows only. The compilation
holds **28** samples called `G-L166(*)`, of which **21 are the hedged `G-L166*` form**. The
`haplotree.info` query above searched the literal string `G-L166` and therefore never matched
them. An earlier reading of this — that haplotree "puts zero Oberbipp samples at L166", so was not
the source of the claim — was **wrong**; the Oberbipp males are present as `G-L166*`. The 7/21
split is exactly the confirmed/provisional split described below.*

**These are third-party assignments and are not evidence until tested at read level.** Two
independent reasons for caution, both observed here rather than assumed:
- Querying the same site for `G-PF3239` returns only **6** samples, of which just `MX182` and
  `MX183` are from Oberbipp. The paper assigns **7** Oberbipp males `PF3239`, and assigns `MX183`
  to `FGC7739/Z6488` — one node *shallower* than haplotree places it. So haplotree both under-covers
  the published set and moves at least one sample deeper than its source.
- Aesch and Muttenz, which carry 10 and 2 published `PF3239` individuals, do not appear at all.

Aggregator coverage is therefore partial and its labels drift downward — the same failure mode
this project documents for the Iceman. Treat the table above as a to-test list, not a result.

**Provenance of these assignments: haplotree.info does not state one.** Checked `index.php` and
the sample listing on 2026-07-25; neither carries a citation, curator, or methodology statement,
and the version string "All Ancient DNA v. 2.07.26" matches no upstream compilation's scheme. The
recoverable trail is the sample IDs: `I#####` are Reich Lab accessions and therefore in AADR;
`MX###`/`SX##`/`TU876` are Furtwängler 2020 / Tübingen; `UNTA58_68Sk1` is Unterer Talweg 58–62,
Augsburg (Lech Valley); `NOE001` and `E09538` follow other schemes and are not yet traced.

The canonical route is the **AADR** (v66.0, 2026-04-13; 19,119 sequences / 17,634 individuals),
whose `.anno` file carries per-sample metadata extracted from the reporting papers, including the
publication. That gives sample → paper → ENA/SRA → reads, which is the same chain used for the
Swiss samples here and the only one that supports re-deriving a call.

**Caution specific to `I5118`:** it also appears on DNAGENICS asserting `G-L166` — the aggregator
already documented above as presenting `MX211` (published `PF3239`) as `G-L166`, crediting
PMID 33135465 (a pathogen target-enrichment methods paper) and describing its whole method as
"We use the latest phylotree for YDNA haplogroup classification and data". Two aggregators agreeing
is not independent confirmation when neither shows its evidence and one may be ingesting the other.
`I5118`'s source publication has **not** been identified; do not infer it from the ID prefix.

**The Horgen Claim Checked Against Furtwängler 2020 Suppl. Tables 1 and 5 (2026-07-25)**
Read from the *corrected* Supplementary Information
(`static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15560-x/MediaObjects/41467_2020_15560_MOESM1_ESM.pdf`),
not from a search summary. Three of the circulating claim's premises fail against its own source.

**1. `PF3239` is not Horgen-specific — Aesch has more of it than Oberbipp.** Individuals assigned
terminal `PF3239` across the whole study, by site:

| site | n | culture label in the table |
|---|---|---|
| Aesch (CH) | **10** | *(none given)* |
| Oberbipp Horgen (CH) | 7 | Horgen |
| Muttenz (CH) | 2 | *(none given)* |
| Niederried Ursisbalm (CH) | 1 | *(none given)* |
| **total** | **20** | |

The claim that this lineage is not dominant in cultures contemporary with the Iceman "except for
the Horgen Culture" is contradicted by the single dataset it rests on. Only Oberbipp carries the
string "Horgen" in the site column at all.

**2. The study's only `L166` is a different site from the seven, and the publication assigns it no
culture.** `L166` occurs exactly once in the entire supplement: `SX10`, `G2a2a1a2a1a`, the
unborn/newborn baby of `SX8`, at **Rapperswil Zürichstrasse**. The seven Oberbipp males are
`PF3239`. So the one genuine ancient `L166` in this dataset is precisely the sample the claim does
not cite.

*Correction (2026-07-25, same day): an earlier version of this note said `SX10` "is not a Horgen
sample". That overstated the evidence.* The paper describes Rapperswil only as a stone cist burial
dated 2695–2481 BC and gives it **no cultural attribution anywhere** — the Horgen label appears in
the site column for Oberbipp and for no other site. `haplotree.info` does attribute `SX10` to the
Horgen culture. Whether Rapperswil is culturally Horgen is an archaeological question this
publication does not answer, and the date sits at the Horgen/Corded Ware transition. The defensible
statement is that the paper does not call it Horgen, not that it isn't. If the attribution holds,
then a genuine Horgen `L166` does exist — one unborn infant, which is still not "7 of 10".

**3. The Oberbipp males are a kin group, so they are not seven independent observations.**
Supplementary Table 5 (PMR / READ / lcMLkin) lists first- and second-degree pairs among them,
including `MX187`–`MX212` (first degree, both in our test set) and `MX187`–`MX209` (second
degree). Male–male first-degree pairs are father–son or full brothers and therefore carry the
same Y chromosome by descent. "7 of 10 carry the clade" is closer to "one patrilineage, sampled
repeatedly"; the paper's own text describes kinship over three generations in the paternal line.

**4. …and that same table shows the published terminal labels are coverage floors.** Table 5
carries a `same Y HG` column whose values include **`same clade`**, used for related pairs whose
published terminal labels differ. Among Oberbipp: `MX150` (`L91`, `G2a2a1a2`) is first-degree with
`MX187` (`PF3239`, `G2a2a1a2a1`); `MX183` (`Z6488`) is first-degree with `MX211` (`PF3239`);
`MX209` (`PF3239`) is first-degree with `MX219` (`PF3147`, `G2a2a`). These men must share a Y
chromosome, so a label difference of one to three nodes between them cannot be phylogenetic — it
is depth. **The Y column of Table 1 must not be read as terminal placement**, which is the same
error this project documents for the Iceman, appearing here inside the source publication's own
data.

The consequence cuts both ways and is worth stating plainly: point 4 is direct published support
for the *premise* of H1 — the `PF3239` stop is coverage-limited, so these men could be `L166`.
The circulating claim could therefore be accidentally correct while every premise it states is
wrong. That is exactly what the read-level test is for. Note also that point 4 applies equally to
Aesch's ten, so "Horgen is special" fails under H1 as well as under H2.

**y_markers_pileup.py Did Not Implement The Registered Call Rules (2026-07-25)**
Found by running the first Swiss sample through it. `y_markers_pileup.py` called any site with
`der > 0, anc == 0` **DERIVED** and any site with `anc > 0, der == 0` **ancestral**, with no depth
threshold and no damage rule. Both contradict the decision rules in
`PREREG_swiss_neolithic_L166.md`, and one of them contradicts the point of writing that document:

> Ancestral call: >=2 independent reads. A site with 1 ancestral read is reported as `low_power`,
> **not** as ancestral.

On MX182 the old logic reported `L166` as **ancestral** off a single read at DP 1. That is H0 — no
power — presented as H2 — tested and negative. The prereg names that specific substitution as "the
single most important commitment in this document", so the tool was violating the registered
protocol on the very first sample it saw.

The rules now live in `site_call()` in `y_markers_pileup.py`, identical to the logic in
`y_sites_pileup.py`, and MX182 correctly reports `low_power_1read_ancestral`.

**Nothing about the Iceman conclusion moves.** Re-running the Iceman BAM through the corrected tool
gives byte-identical output for `L166_defining` (9 DERIVED), `Z6494_exclusion` (3 ancestral) and
`backbone_control` (10 DERIVED) — every call there rests on DP>=2 with unanimous support, so the
stricter rules cannot touch it. Across all committed Iceman tables only **3 calls** change, all
DP-1 deamination-prone transitions and all in the conservative direction:

| file | marker | change |
|---|---|---|
| `iceman_y_ftdna_block_evidence.tsv` | `FGC2274` | DERIVED -> `nocall_damage_prone_1read` |
| `iceman_y_L166_subtree_exhaustive.tsv` | `FT172194` | DERIVED -> `nocall_damage_prone_1read` |
| `iceman_y_L166_subtree_exhaustive.tsv` | `Y126330` | DERIVED -> `nocall_damage_prone_1read` |

`iceman_y_L166_evidence.tsv` and `iceman_y_ftdna_project_terminals.tsv` are unchanged.

All three had already been dismissed as artifacts in prose above ("0 credible derived", "DERIVED
(weak)"), so the analyst's reading was right and only the `call` column was out of step. The useful
consequence is that the L166 subtree scan now tallies **88 ancestral / 1 mixed / 2 damage-nocall /
2 no-coverage — zero derived** as data, rather than showing two DERIVED rows that had to be argued
away in a footnote.

**Ad-hoc MAPQ Table Reproduced By Committed Tool (2026-07-25) — two errors found in the table**
The MAPQ audit below was originally produced by a process that was never committed, which meant the
single most consequential result in the project could not be regenerated by anyone, including its
author. `annotate/y_sites_pileup.py` now does it. Re-run over all 21 sites
(`markers/iceman_novel_candidates_all21.tsv`) against the same Iceman BAM it reproduces
**21/21 verdicts**, REJECT/MARGINAL classification included. Output:
`results/iceman_y_novel_candidates_regen.tsv`.

Five cells differed, and in every case the committed tool is right and the hand-built table wrong:

- **`chrY:11554430` `fwd` was 6, should be 7.** The published row has `fwd + rev = 11` against
  `dp = 12` — its own columns did not add up. The base string is `.AaAAAaAaA^Ja^Na`: one forward
  reference read plus six forward `A`, five reverse. The ad-hoc pass appears to have dropped the
  ancestral read from the forward count.
- **Four `pct_mq0` cells differ by exactly 1%** (`11401880`/`11401881` 43→44, `11554430` 83→84,
  `11578588` 77→78). Raw values are 43.75 / 83.82 / 77.55%. The old pass truncated, the tool rounds.
  No verdict moves; all four are rejections either way.

The reproduction also surfaced something the old table could not represent: **`chrY:11554430` is
mixed**, 1 ancestral G against 11 derived A. That table has no `n_anc` column, so a mixed site there
was invisible by construction. The verdict is unchanged — it is rejected at 84% MQ0 regardless — but
this is exactly the failure mode a hand-built summary hides, and it is why `n_anc`/`n_der`/`n_other`
are now emitted separately from the call.

Two rules from the prereg existed only in prose and are now in the code: the `pct_mq0 >= 30%`
rejection threshold (`site_qc` column) and no-go region flagging (`region` column). The no-go flag is
advisory and deliberately does not reject, because applying it retroactively would move
`chrY:11414525` and `chrY:11667647` — both inside the 11.1-11.7 Mb window — out of the usable set
that this note already published.

**Novel-Branch Candidates Re-derived With MAPQ Filtering (2026-07-25) — supersedes earlier list**
The previous candidate list filtered on depth only. Adding mapping quality changes it substantially:
**10 of 21 candidates fail at 39-88% MAPQ-0 reads**, including two of the three transversions that had
been presented as the strongest damage-immune evidence (`chrY:11176570` at 70% MQ0, `chrY:11578588` at
77%). `chrY:20327833`, added earlier the same day as "marginal", is also rejected at 44%.

Depth alone cannot catch this: these sites have unremarkable DP (4-12) and only reveal themselves in the
MAPQ distribution. The chrY:11.1-11.7 Mb window is particularly bad and should be treated as a
no-go zone alongside 56.69-56.88 Mb and ~26.6 Mb.

Surviving set in `results/iceman_y_novel_branch_candidates.tsv` — **8 usable, 3 marginal**:
- Best single site: `chrY:7885869 T>A`, a transversion at DP 9, both strands, **0% MAPQ-0 and 10/10
  reads at MQ60**. This is the strongest private-variant candidate in the dataset.
- Also clean: `7899558 A>G`, `8808031 C>T`, `10964462 T>C`, `11667647 T>C`, `19647870 C>T`,
  `11414525 T>C` (7% MQ0), `10768171 G>A` (12% MQ0).
- Marginal: `10990649` and `10996925` (no MQ60 reads at all), `13414189 C>G` (transversion, clean
  mapping, but 5/5 reads reverse-strand).

Only one high-confidence transversion now survives, not three. Columns record `pct_mq0`, `n_mq60`, and
forward/reverse read counts so this is auditable rather than a bare verdict.

`CGG017683` was re-tested at the surviving sites and covers 3 of them (`7899558`, `11414525`,
`19647870`) — **ancestral at all three**. The negative result from the earlier comparison holds and is
now based on 5 informative sites rather than 2.

**Remapping Note (correction)**
Full FASTQ is available for these samples (`ERR14752008.fastq.gz` etc.), so a proper remap to hg38 through
`mapping/revert-bam.sh` is possible with **no** reference bias — the bias objection only applies to
re-aligning an already-Y-filtered slice. Remapping is therefore legitimate here, and would additionally
recover reads lost to the depositors' `MQ>=25` pre-filter. It just does not help `CGG017683`, where the
limit is 0.192x depth rather than alignment quality.

**Iceman Novel-Branch Candidate SNPs (superseded — see the MAPQ-filtered re-derivation below)**
Kept for the reasoning; the numbers here are stale.
- First pass (VCF panel, depth filter only): 17 novel sites, 10 rejected as collapsed-repeat pileups at
  `chrY:56.69-56.88 Mb` (Yq heterochromatin / PAR2 boundary, DP 28-269 vs a ~8x chrY mean),
  leaving 7 usable of which 3 were called transversions.
- Both figures changed on re-derivation: the GFF3 index gives a different novel set than the VCF panel,
  and adding a MAPQ filter rejects two of those three transversions. Current numbers live in
  "Novel-Branch Candidates Re-derived With MAPQ Filtering".
- The intent stands: this is the target list for defining a subclade below `G-L166` if a branch-mate
  ever appears.

**Marker-Level Artifacts in Repetitive chrY (2026-07-25)**
- `Z6519` (chrY:21049898) scores derived in the callset but is an artifact:
  at `-q 0` the site is 9 reads, 7 ancestral — and all 7 ancestral reads are MQ 0.
  Only the 2 derived reads survive a MAPQ filter, so the site flips with the filter setting.
- Every other marker sharing its label (`G2a2a1a2a1a1b2~`, ~24 markers) is ancestral at DP 4-13.
- Check MQ distribution at any single-marker "deep branch" hit before believing it.

**VCF vs gVCF Path-Rank Note**
- gVCF can slightly change path ranking compared with VCF because explicit non-variant blocks reduce `nocall` counts at marker positions.
- Low-score candidates present only in gVCF top-path output should be treated as weak evidence unless they also gain additional derived-marker support.
