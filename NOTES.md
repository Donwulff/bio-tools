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
  P(any read at a given site) = 4%. All reads pre-filtered to **MQ>=30** (read from the deposited
  BAM header 2026-07-26; this line previously said MQ>=25, which was wrong).
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
- Aesch and Muttenz, which carry 10 and 2 published `PF3239` individuals, ~~do not appear at all~~
  *do not appear in that query*.

*Correction (2026-07-26): the struck clause was wrong, and wrong in the same way as the `G-L166`
literal-string search corrected immediately above. Aesch and Muttenz are present in the compilation
in quantity — as `Aes1`–`Aes25` and `RA42`–`RA64`, IDs that no search for the strings "Aesch" or
"Muttenz" will match. Between them they carry **15 `G-L166*` rows**, which are the Test C cohort
registered in `PREREG_testC_aesch_muttenz.md`. What is true is the narrower statement: neither site
appears under a `G-PF3239` query, because the compilation does not label them `G-PF3239`.*

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

**Remapping Note (correction)** — *itself corrected 2026-07-26, see the strike below.*
Full FASTQ is available for these samples (`ERR14752008.fastq.gz` etc.), so a proper remap to hg38 through
`mapping/revert-bam.sh` is possible with **no** reference bias — the bias objection only applies to
re-aligning an already-Y-filtered slice. Remapping is therefore legitimate here, and ~~would additionally
recover reads lost to the depositors' MQ pre-filter (`MQ>=30` for `CGG017683`)~~. It just does not help `CGG017683`, where the
limit is 0.192x depth rather than alignment quality.

**The struck clause is wrong and contradicts "Every Sample In This Project Came From A Mapped-Only
Deposit" below, which is authoritative.** Re-mapping cannot recover reads lost to the depositor's MQ
pre-filter, because the ENA FASTQ is *generated from the already-filtered BAM*: for `ERR14752008` the
FASTQ `read_count` is 27,290,045 and the BAM's mapped total is 27,290,045 with 0 unmapped. The FASTQ
**is** the filtered BAM. Every read the depositor's `samtools view -q 30` discarded is absent from
both. What re-mapping does buy is real but narrower: native hg38 placement instead of an
hg38→hg19 liftover of marker coordinates, our reference, and our own MAPQ and damage rules applied
to the reads that survived. The conclusion — that re-mapping does not rescue `CGG017683` — is
unchanged, but the reason is depth *and* an irrecoverable read set, not depth alone.

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

**SWISS NEOLITHIC L166 TEST — RESULT (2026-07-26)**
All 9 samples mapped (`mapped=9 skipped=0 failed=0`), genotyped by
`annotate/y_genotype_batch.sh`, verdicts by `annotate/y_prereg_verdict.py`.
Tables: `results/swiss/`.

| sample | chrY DoC | verdict |
|---|---|---|
| MX210 | 0.1434x | `CONFLICT` |
| MX211 | 0.0769x | `H0_no_power` |
| MX212 | 0.0561x | **`H2_ancestral`** |
| MX187 | 0.0545x | `H0_no_power` |
| SX10 | 0.0485x | `H1_derived` |
| MX213 | 0.0416x | **`H2_ancestral`** |
| MX182 | 0.0389x | `H0_no_power` |
| MX209 | 0.0295x | `H0_no_power` |
| SX8 (female control) | 0.0009x | `H0_no_power` |

`SX8` at 0.0009x against 0.0295-0.1434x for the males is a clean negative control.

**Only 3 of the 9 registered markers are reachable in this dataset.** `FGC5696`, `Z6208`,
`S19530`, `FGC5721` and `Z6516` have **zero reads in all eight males** — systematic, not
coverage variation. `Z6287` appears once. The registered 9-marker primary test is a 3-marker
test in practice, and that was not foreseeable from depth.

**The three reachable markers disagree systematically:**

| marker | class | reads | derived | ancestral |
|---|---|---|---|---|
| `Z6219` | C>T transition | 10 (6 samples) | **10** | 0 |
| `L166` | C>A transversion | 11 (6 samples) | 0 | **11** |
| `L167` | T>A transversion | 5 (4 samples) | 0 | **5** |

`Z6219` is derived in every sample covering it; `L166`/`L167` are ancestral in every sample
covering them, *in the same samples*. `Z6219` is therefore **not diagnostic for `L166`** — it is
derived at some node above it, or recurrent. This is a post-hoc observation, not a registered
rule, and the verdicts above are reported as the registered rules produced them.

Consequences:
- `MX210`'s `CONFLICT` is entirely `Z6219` versus the two transversions. On the damage-immune
  markers alone it is ancestral, 3/3 and 2/2.
- **`SX10`'s `H1_derived` rests on `Z6219` alone** (DP 2, one read terminal), with *zero* reads at
  `L166` and `L167`. The positive control therefore does not confirm that this panel detects
  `L166` where it exists; it confirms that `Z6219` is derived in a sample published as `L166` —
  and `Z6219` is equally derived in samples demonstrably ancestral at `L166`.

**Damage cannot explain the ancestral calls.** Deamination runs C->T and G->A, so it cannot
manufacture an ancestral `C` at the C>A site `L166`, nor an ancestral `T` at the T>A site `L167`.
All 16 transversion reads are at 0% MQ0.

**Bottom line: 16 damage-immune reads across 6 Oberbipp males, every one ancestral at the L166
locus, zero derived reads anywhere.** Two samples reach the registered `H2_ancestral` threshold
and a third (`MX210`) does so on the transversions alone. Under the registered rules this does
not license a population-level statement — these men are a kin group, `Family A` covers 9 of
them, and no proportion may be computed over them. What it does say is that the compilation's
`G-L166*` relabelling is not supported at read level in any sample where the defining
transversions are observed, and that the compiler's own `G2a2a1a2a1(a)` hedge was correct.


**Extension To 15 Samples: Four Independent Lineages, Still Zero Derived Reads (2026-07-26)**
Six further Oberbipp runs staged from `PRJNA608699` and genotyped under `PREREG` Amendment 2.
Results in `results/swiss15/`; the 9-sample tables in `results/swiss/` are left untouched as the
record of what was published from the first pass.

*The kin structure is richer than the earlier note assumed.* `Kinship-Notes` in the compilation
names **two** families at Oberbipp, not one: `MX213` is `Family B (3 members)`, not merely
"outside Family A". `MX204`, `MX210` and `MX299` carry `n/a`. So the assemblage yields **four**
Y-chromosome observations that are independent of each other:

| lineage | evidence at `L166` | at `L167` |
|---|---|---|
| **Family A**, 9 members pooled | **ancestral 7/7 reads** (5 contributors) | **ancestral 3/3** (3) |
| **Family B** (`MX213`) | ancestral 3/3 | ancestral 1/1 |
| `MX210` (`n/a`) | ancestral 3/3 | ancestral 2/2 |
| `MX204` (`n/a`) | ancestral 1/1 | ancestral 1/1 |

**Across all 15 samples: 21 reads at the two `L166`-defining transversions, every one ancestral,
zero derived.** All at 0% MQ0, all `site_qc pass`. Both sites are transversions, so deamination
cannot manufacture the ancestral allele at either. Per-sample verdicts add `MX183` to the
`H2_ancestral` list (`MX183`, `MX212`, `MX213`); `MX210` stays `CONFLICT`; no Oberbipp sample
reaches `H1` at any depth.

**The Family A pool behaved exactly as prediction E registered.** It sharpened the conflict rather
than resolving it: `Z6219` derived 3/3, `L166` ancestral 7/7, and the five markers with zero
coverage in every individual (`FGC5696`, `Z6208`, `S19530`, `FGC5721`, `Z6516`) stayed at zero
even pooled. The pool's value is the positive control, which is now unambiguous for the patriline
— `P15`, `M3308`, `PF3147`, `L91` (12 of 13 reads, transversion), `PF3239`, `Z6043`, `P287` all
`DERIVED`, no ancestral read anywhere on the backbone. The published chain is recovered; the node
below it is not. **This pool is one observation about one chromosome and may not be written as a
proportion over the Oberbipp men.**

**The `Z6219` localisation test returned H0, as pre-registered.** `MX204` was the one sample in
the study that could have settled it — published `G2a2a1a2a` (`Z6488`), no kinship record, 0.0320x
chrY. It has no reads at `PF3239`, so it fails the registered admissibility gate and is **not** an
outgroup; it also has no reads at `Z6219`, so the question could not have been answered even had
it qualified. `MX299` is empty throughout. **Whether `Z6219` is a mis-levelled catalogue entry or
a genuine node between `Z6488` and `L166` remains open, and this dataset cannot close it.** The
outgroup has to be a G2a male demonstrably ancestral at `PF3239`, and the Oberbipp assemblage does
not contain one with usable chrY depth.

*What `MX204` did contribute* was not what it was staged for: a fourth independent lineage,
ancestral at both defining transversions. One read each is `low_power_1read_ancestral` under the
standing rules — `H0`, not `H2`, and it is not counted as evidence — but it is recorded because
the rule that excludes it is the same rule that would have excluded a convenient result.

**`MX203` is almost certainly female.** 0.0001x on chrY, below the known-female control `SX8` at
0.0009x, against 0.142 autosomal coverage — so the library is fine and the chrY is absent. It is
the ninth `Family A` member by the published kinship note and is listed in
`markers/family_a_members.txt` for completeness; it contributes no reads and is dropped by the
pooling script on its own merits rather than by special-casing.

**No shared novel variant.** Of the six added samples only `MX183` has any read at the 8-site
Iceman novel panel — one read at `chrY:7885869`, and it is **ancestral**. The panel remains at
zero informative sites, as registered.

**Two documentation defects found and fixed while staging.** `fetch_ena_runs.sh` wrote
`manifest.tsv` with `>`, so staging a second subset into the same directory would have erased the
provenance record for the nine files already there; it now merges on `run_accession`. And the
`RUNLOG` staging command named `PRJEB36959`, which ENA resolves to *"Partial sequence of 16S
Leptospira borgpetersenii and Leptospira interrogans"* and which returns zero read runs — the
study is `PRJNA608699`. It would have failed loudly rather than fetched wrong data, so no result
is affected, but it did not reproduce. Related: the earlier session queried
`haplogroup.info/all-ancient-dna.txt` without keeping a copy, so the claims resting on it were not
re-runnable; a dated copy is now kept at `/mnt/AncientDNA/all-ancient-dna.2026-07-26.txt`.

---

## The Unhedged `G-L166` Set: The Label Fails In Both Directions (2026-07-26)

Pre-registered in `PREREG_unhedged_L166.md`, committed `e37b8ab` before any read was staged.
Tables in `results/unhedged/`. Five individuals, three studies, three countries.

**The set is six people, not seven.** `UNTA58_68Sk1` (Mittnik 2019) and `E09538` (Olalde 2018)
share radiocarbon lab number **MAMS-29075** (3870±30 BP), mtDNA `J1c` and coordinates, and the
Mittnik identifier decodes to the Olalde colloquial description — "Feature 68 Skeleton 1" is site
UNTA58, feature 68, skeleton 1. One man, two aggregation IDs. `SX10` was already done, leaving
five to stage.

**The unhedged label does not verify, and the failures are deep, not marginal.**

| sample | site | chrY DoC | `L166` | `L167` | `FGC5696` | verdict |
|---|---|---|---|---|---|---|
| `I14677` | Serra Crabiles, Sardinia | 0.2973x | **anc 11/11** | anc 5/5 | anc 5/5 | `H2_ancestral` |
| `I14678` | Serra Crabiles, Sardinia | 0.3514x | **anc 8/8** | anc 5/5 | anc 5/5 | `H2_ancestral` |
| `I15942` | Anghelu Ruju, Sardinia | 0.1013x | anc 2/2 | anc 1 (low power) | — | `H2_ancestral` |
| `I5118` | Mezőcsát-Hörcsögös, HU | 0.1406x | **der 1** | — | **der 1** | `H1_derived` |
| `UNTA58_68Sk1` | Augsburg, Bavaria | 0.0597x | — | — | — | `H0_no_power` |

All calls at 0% MQ0 and `site_qc pass`. **42 ancestral reads across the three L166-defining
transversions in the three Sardinians, zero derived.** These are transversions, so deamination
cannot manufacture the ancestral allele. This is not a coverage floor: `I14677` carries 11 reads
at `L166` itself, deeper than any Oberbipp sample managed anywhere.

So the outcome is **H3, the split case registered in advance** — and the split does not run the
way the hedged/unhedged distinction predicts. Two of three Sardinian individuals carry the
*unhedged* label at depth sufficient to refute it outright.

**Counting independent lineages, per the registered rule.** `I14677` and `I14678` are from the
same tomb (Serra Crabiles, Tomb III Cella A, identical coordinates) and were registered in advance
as **one** observation unless their reads showed otherwise. Their genotypes are identical at every
informative marker, which is consistent with one patriline and does not distinguish it from two
men of the same clade — so the rule stands and they count as one. That gives **two independent
lineages falsifying the label** (`I14677`/`I14678`, `I15942`) against **one verifying it**
(`I5118`), with one untestable.

**`I5118` verifies, and it is the interesting one.** Hungary, 3300–3000 BCE, near-contemporary
with the Iceman. Derived at **two** independent L166-defining transversions — `L166` and
`FGC5696` — one read each, with no conflicting read, which is a `DERIVED_1read_transversion` call
under the standing rules. Backbone `P15`, `PF3147`, `L91`, `PF3239`, `Z6043` all DERIVED with no
ancestral read. It is the only sample here that is `L166`-derived, and it is the only
non-Sardinian, non-Bell-Beaker individual in the set. **Its depth is 1 read per site.** The rule
permitting that was fixed before any of this data existed and is not being revisited, but the call
is a single molecule at each site and is reported as exactly that.

**`UNTA58_68Sk1` is H0 and its library is not what the metadata says.** Zero reads at all nine
L166-defining markers. It retains 14,652 of 15,092 chrY reads at MAPQ≥25 (97%), against 41–43% for
the four capture libraries — the signature of shotgun rather than capture, corroborating the
contradictory `Targeted-Capture`/`RANDOM` deposit annotation flagged in the prereg. Several of its
sites are `MARGINAL_no_unique_reads`, i.e. what little is there is duplicate.

**No `L166*` call is possible here, as registered.** `Z6494` is off-panel; only `Z6215` (a
transition) carried reads. `I14678` is ancestral 2/2 there, so it is not `G-Z6494`, but the
exclusion cannot be made to the standard applied to the Iceman.

**No shared novel variant.** Nothing on the 8-site Iceman panel reaches an informative call in any
sample. `I14678` has 2 mixed reads at `chrY:11414525`, which is inside the declared 11.1–11.7 Mb
no-go window and flagged `MARGINAL_single_strand`; `I15942` has a `REJECT_mapq(100%_MQ0)` at
`10964462`. Both are non-findings, as predicted.

### The `Z6219` conflict now has an answer — reached post-hoc, and labelled as such

`PREREG_swiss_neolithic_L166.md` Amendment 2 §C left two explanations for `Z6219` being derived in
every Oberbipp male while `L166` was ancestral in the same men: **(i)** the catalogue mis-levels
it and it sits at or above `PF3239`; **(ii)** it marks a real node between `PF3239` and `L166`.
That test returned H0 exactly as its registered prior predicted.

The Sardinians settle it:

| | `PF3239` | `Z6219` | `L166` |
|---|---|---|---|
| `I14677` | **DERIVED 5/5** | **ancestral 4/4** | ancestral 11/11 |
| `I14678` | **DERIVED 3/3** | **ancestral 2/2** | ancestral 8/8 |
| Oberbipp (15 samples) | DERIVED | DERIVED 10/10 reads | ancestral 11/11 |
| — of which `MX210` alone | DERIVED | **DERIVED 4/4** | **ancestral 3/3** (`L167` anc 2/2) |
| Iceman | DERIVED | DERIVED | DERIVED |

*Correction to the Oberbipp row (2026-07-26).* `10/10` is a correct **read** tally but overstates the
**call**-level support, and the distinction matters because `Z6219` is `C>T` — the one direction
deamination manufactures. Broken out: `MX210` 4/4 derived (a `DERIVED` call), `SX10` 2/2 derived
(`DERIVED`, but no `L166` coverage, so it constrains nothing here), and `MX187`/`MX211`/`MX212`/
`MX213` at one derived read each — every one of which the registered rules score
`nocall_damage_prone_1read`. Three of those four are the Family A patriline and pool to a single
observation. **The "above `L166`" constraint therefore rests principally on `MX210`**, the only individual
carrying both a clean derived call at `Z6219` and clean ancestral calls at two `L166`-defining
transversions, and which is outside Family A.

*Superseded 2026-07-26 by the read-terminus check below: the four single reads are not damage
either, and the constraint is considerably stronger than this paragraph allows.*

An individual derived at `PF3239` and ancestral at `Z6219` places `Z6219` **below** `PF3239`,
which **refutes (i)**. Oberbipp places it above `L166`. The consistent topology is
`PF3239 → Z6219 → L166`, with the Sardinian branch splitting off before `Z6219` and the Oberbipp
patriline between `Z6219` and `L166` — i.e. **explanation (ii), a real node**, and the Oberbipp
men sit one node below their published `PF3239` call, closer to the Iceman than the publication
states and still not `L166`.

The direction of the damage argument favours this reading. `Z6219` is C>T, so deamination converts
ancestral to derived; it **cannot** manufacture the ancestral `C` observed 4/4 and 2/2 here. The
weaker link is `PF3239`, also C>T, where the derived call could in principle be damage — but both
samples are independently `L91`-derived on transversions (6/6 and 10/10), so they are below `L91`
regardless, and `PF3239` at 5/5 and 3/3 with 0% MQ0 is not a plausible damage artefact.

**This inference was not pre-registered, and the prereg for this very dataset says these samples
could not do it.** That claim was wrong, and it was wrong for an instructive reason: it assumed
all five samples were genuinely `L166`-derived *because that is what their published labels say* —
the very thing under test. The failure of the label is precisely what produced the outgroup I had
argued did not exist. Registering a prediction on the basis of the labels, in a study whose
purpose is to check the labels, is a small circularity worth naming.

The registered admissibility gate (Amendment 2 §C) required an outgroup **ancestral** at `PF3239`,
and discriminated in the opposite direction. The configuration actually observed — **derived** at
`PF3239`, ancestral at `Z6219` — is equally valid logically but is not the rule that was written
down. Per the standing caveat in that amendment, this is therefore **a hypothesis about tree
topology, well supported but resting on one marker, and it requires a pre-registered confirmation
before it is treated as settled.** It is not being written into `markers/L166_defining.txt` on the
strength of a post-hoc reading.

### What this does and does not change

- The circulating "7 of 10 Oberbipp are `G-L166`" claim was already refuted by 21 ancestral reads
  across four independent lineages. Nothing here weakens that.
- The hedged/unhedged distinction in the compilation **does not carry the read-level meaning this
  project inferred for it**. Unhedged `G2a2a1a2a1a` samples are ancestral at `L166` at depths of
  8 and 11 reads. The inference is not merely unconfirmed, it is contradicted, and the
  reading of the parenthetical convention in earlier notes should be treated as unsupported.
- Genuine `G-L166` individuals other than the Iceman **do exist** — `I5118`, in Hungary, at
  near-contemporary date — but on one read per site at two markers.
- Nothing here speaks to where anyone came from. The geographic spread invites that reading and
  the prereg refuses it in advance.

### The compilation's own evidence column already disagreed with its label column (2026-07-26)

Asked afterwards where these labels came from, and whether the answer is motivated reasoning
somewhere upstream. Checked `all-ancient-dna.2026-07-26.txt` directly rather than reasoning about
it. The file carries, per row, a `SNP-positive` list of every derived call with the clade each SNP
defines. Comparing that field against the `ISOGG2019` label field, for these six:

| sample | deepest node in `SNP-positive` | `ISOGG2019` label | `Responsible-SNP` | read-level verdict |
|---|---|---|---|---|
| `I14677` | `G2a2a1a2a1` (PF3239, FGC5666, Z6130, +2) | `G2a2a1a2a1a` | Kolgeh | ancestral 11/11 at `L166` |
| `I14678` | `G2a2a1a2a1` (PF3239, Z6130, Z6277, +2) | `G2a2a1a2a1a` | Kolgeh | ancestral 8/8 |
| `I15942` | `G2a2a1a2a1` (PF3239, Z6130, Z6277, +2) | `G2a2a1a2a1a` | Kolgeh | ancestral 2/2 |
| `I5118` | `G2a2a1a2a1a` — **incl. `L166:23989884C->A`** | `G2a2a1a2a1a` | *(blank)* | **derived** at `L166` |
| `E09538` | `G2a2a1a2a1a` (FGC5671^, Z6219 — *not* `L166`) | `G2a2a1a2a1a` | *(blank)* | — |
| `UNTA58_68Sk1` | **empty — no positive SNPs at all** | `G2a2a1a2a1a` | Kolgeh | H0, no power |

**The three samples the reads falsified are exactly the three whose own evidence field never
reaches `L166`.** It stops at `PF3239` in all three. The one sample that verified is the one whose
evidence field lists `L166` itself, at the correct coordinate and as the correct transversion
(`C->A`) — the same call this project made independently from reads. The concordance between the
compilation's evidence column and the pileups is complete; the disagreement is entirely between the
compilation's evidence column and its own label column.

So no suppressed or unexamined data is required to explain this. Someone did the work, recorded it,
and the label field was then filled one node deeper than the recorded evidence supports. That is the
**same structural failure** as the Horgen hedge, which survives in `ISOGG2019` and is lost one column
over in `Y-Haplotree-Variant`, and the same as the Iceman's own `G-Z6208` — a derived SNP read as a
terminal placement. Three occurrences, three different columns, one mechanism: **the qualifier and
the evidence live in a different field from the label, and only the label propagates.**

`UNTA58_68Sk1` deserves separate mention: it carries the unhedged label with an entirely empty
`SNP-positive` field, and its `Assessment` is `QUESTIONABLE (Xcontam=[0.013,0.029])`. Its duplicate
row `E09538` is where the Y evidence sits. Neither row's evidence includes `L166` itself.

**Source publications, and what they do not tell us.** `Source` gives `FernandesNatEcolEvol2020`
for the three Sardinians and `OlaldeNature2018` for `I5118` and `E09538` — both Reich-lab papers.
This does **not** make the labels theirs. `haplogroup.info` states its own method (Yleaf v2 plus
hobbyist BAM analyses) and credits `Responsible-SNP` per row; `Source` names the paper the *sample*
came from, not the origin of the Y call. What Fernandes 2020 and Olalde 2018 themselves assign to
these individuals has **not been read**, and AADR's `.anno` Y-column resolution remains unchecked
(no copy on this machine). Until one of those is done, nothing here can be attributed to any
publication or lab, and it should not be.

Note also which direction the error runs in the cases that *have* been checked: Furtwängler 2020
assigns the Oberbipp males `PF3239`, and the aggregator moved them deeper. There the publication was
the conservative layer. The Iceman's `G-Z6208` shows publications are not immune — but in this
corpus the over-calling has been downstream of them every time it has been traced.

**Corollary for `Z6219`.** The file assigns `Z6219` to `G2a2a1a2a1a`, i.e. the `L166` node itself.
The Oberbipp series is `Z6219`-derived 10/10 and `L166`-ancestral 11/11, which contradicts that
placement on its own, independently of the Sardinian data above. If `Z6219` sits between `PF3239`
and `L166` as the reads indicate, then `E09538`'s label rests on a marker that never reaches `L166`
either. This remains post-hoc and still requires its own pre-registration.

### Tree-version drift is real, and the compilation flags it with an undocumented caret (2026-07-26)

Raised as a hypothesis: the labels may not be errors at all, but fossils of an older tree in which
`PF3239` sat inside `G-L166`. Checked against `resources/marker_index.tsv.gz` and the compilation.

**Version drift is demonstrably present inside a single catalogue file.** Markers at the same node
carry two different ISOGG longhands:

| longhand | markers | yfull_node |
|---|---|---|
| `G2a2a1a2a1a` | `L166`, `L167`, `Z6219`, `Z6287`, `Z6516` | `G-L166` |
| `G2a2a1b1a` (superseded numbering) | `Z6134`, `FGC5671`, `FGC5696`, `FGC5721` | `G-L166` |

Both name the `L166` node. `Z6134` and `Z6516` additionally carry the free-text comment
`Equiv. to L166`.

**The `PF3239`/`L166` boundary has genuinely moved.** The compilation assigns `FGC8009` and `PF3178`
to `G2a2a1a2a1a` (`L166`); the catalogue assigns both to `G-PF3239`. Two markers, opposite sides of
the boundary, two sources. So the hypothesis is not idle.

**The caret means something.** The compilation suffixes some SNP names with `^`. No legend has been
located, but across the 15 markers checked here the correlation with the catalogue is exact:

    ^ present  : FGC2271, PF3244, PF3247, S11565, PF3178, FGC8009 (isogg "not listed"/"unknown")
                 FGC5696, FGC5671, Z6134                          (isogg = superseded G2a2a1b1a)
    ^ absent   : PF3239, Z6130, Z6277, FGC5666                    (isogg = current G2a2a1a2a1)
                 L166, L167, Z6219                                (isogg = current G2a2a1a2a1a)

`^` ⟺ *the marker has no current ISOGG longhand*. 15/15. This is an **inference from correlation**,
like the asterisk reading above — do not present it as documented. But if it holds, the compilation
is flagging version-unstable placements per SNP, and the flag is again in a field nobody reads.

**Applying the caret to the six.** Restricting to unmarked (current-longhand) evidence:

- `I14677`, `I14678`, `I15942` — deepest unmarked support is `PF3239`/`Z6130`/`Z6277`/`FGC5666`,
  all `G2a2a1a2a1`. **No unmarked `L166`-level SNP at all.** The reads agree.
- `I5118` — unmarked `L166`-level support is `L166` itself and `Z6219`. The reads independently
  confirm `L166` (and `FGC5696`). Stands either way.
- `E09538` — unmarked `L166`-level support is `Z6219` **only**; `FGC5671` is carets-and-superseded.
  If `Z6219` sits above `L166`, this individual's entire claim to the node evaporates.

**Does drift explain the three failures?** Partly, and better than the "label filled deeper than the
evidence" reading offered above — but not on its own. Against it: the compilation's *own row* uses
`G2a2a1a2a1` for `PF3239` and `G2a2a1a2a1a` for `L166`, distinctly, in the same line. Within the
file as it stands they are not equivalent. For it: the label field and the evidence field need not
have been generated at the same time from the same tree, and a stale label beside a refreshed
evidence list would produce exactly this. **Not resolved.** Testing it needs a dated ISOGG tree
series (2016–2019), which is not held locally.

### Why `Z6219` could not have been found in a modern-only tree

The catalogue lists `Z6219` as `isogg G2a2a1a2a1a` / `yfull G-L166` — i.e. inside the `L166`
equivalence block. Splitting an equivalence block requires a sample **on the intermediate branch**:
derived at `Z6219`, ancestral at `L166`. A tree built from living testers can only contain such a
sample if that lineage still has living male-line descendants. If it does not, the block stays
merged permanently, no matter how many moderns are sequenced.

This is a structural capability of ancient DNA rather than a defect in the modern trees, and it
generalises: **every equivalence block is a hypothesis that no intermediate lineage survived.**
Ancients are the only way to test it.

Caveat against over-reading the extinction angle: Y-lineage extinction is the *default* outcome of a
neutral branching process, not a signal of anything. Most patrilines leave no male-line descendants
within a few generations. "Absent from moderns" is therefore the expected state of an arbitrary
Neolithic lineage and is **not** evidence of demographic replacement, violence, or anything else.
Note also that `G-L166` itself is *not* extinct — FTDNA has L166 testers
(`results/iceman_y_ftdna_project_terminals.tsv`). Whatever is missing is sub-branch structure, not
the clade.

## Every Sample In This Project Came From A Mapped-Only Deposit (2026-07-26)

Raised as a suspicion — that ancient deposits are probably filtered to mapped reads, since the
libraries are mostly incidental metagenome. Checked, and it is true of **19 of our 20 own-mapped
samples**. This qualifies the "we did our own mapping so they are comparable" claim in a way that
should have been established before the first BAM was made.

**Every ENA "FASTQ" used here is ENA-generated from a submitted BAM.** `submitted_format=BAM` for
`ERR14752008` (Crimea), `ERR2207344` (Olalde), `ERR3518170` (Mittnik) and `ERR3800863`
(Furtwängler) — the four checked, spanning all four studies. For `CGG017683` the arithmetic is
exact: ENA's fastq `read_count` is 27,290,045 and the BAM's mapped total is 27,290,045, unmapped 0.
The fastq *is* the filtered BAM.

**Confirmed empirically against our own output** (`samtools idxstats`, unmapped fraction of the
BAMs this repo produced):

| sample set | unmapped fraction |
|---|---|
| all 15 Swiss (`MX*`, `SX*`) | 0.000 (2–138 reads each) |
| `I14677`, `I14678`, `I15942`, `I5118` | 0.000 |
| **`UNTA58_68Sk1`** | **0.434** |

Nineteen inputs contain essentially nothing that fails to map, which is only possible if the
depositor had already discarded everything that did not map. `UNTA58_68Sk1` is the sole exception,
and this is now the **third** independent signal that its deposit is a different kind of object from
the others: 97% of chrY reads at MAPQ>=25 against 41–43% for the capture libraries; a contradictory
`Targeted-Capture`/`RANDOM` annotation; and now 43% of the library not human-mappable at all, which
is what an unfiltered shotgun aDNA library actually looks like.

**What our re-mapping does and does not achieve.** It is genuine and worth having: read placement,
MAPQ, damage-aware calling and the reference build are ours, applied identically across samples. But
the *read set* is not ours. It was selected by each depositor's aligner against each depositor's
reference. **We can re-place any read some earlier mapping already accepted; we can never recover a
read that earlier mapping rejected.** Where the depositor used GRCh37 this is a build-specific
ascertainment sitting inside the data — GRCh38 substantially revised chrY, so sequence that exists
only in GRCh38 cannot be represented no matter what we map to. The bias has a known direction (loss
only) but an unknown magnitude, and it cannot be corrected from the deposited data.

**The one depositor pipeline that is publicly readable.** `CGG017683`'s BAM header, in order:

    AdapterRemoval (reads are *.fastq.truncated.gz)
    bwa samse            -> hs37d5
    samtools view -F 4 -q 30
    samtools view -F 1024 -q 30 -b ... 1..22 X Y MT     (drops hs37d5 decoys/unplaced)
    picard MarkDuplicates REMOVE_DUPLICATES=true
    GATK IndelRealigner 3.7
    samtools calmd -Erb  -> hs37d5

Three things follow. The MQ cut is **30, not 25** — the earlier note in this file was wrong and is
corrected above. The `-q 30` was applied **in GRCh37 coordinates**, so any site multi-mapping in
GRCh37 but unique in GRCh38 has already been stripped. And `calmd -Erb` means **extended BAQ was
computed**; whether that capped `QUAL` in place or was written as a `BQ` tag depends on `-A`, which
is absent, so the effect on our `--min-bq 20` threshold is **unverified and needs checking before
this sample's base qualities are compared with any other sample's**.

The submitted BAMs for the Swiss, Sardinian, Hungarian and Bavarian runs are **not** distributed
(`submitted_ftp` is empty despite `submitted_format=BAM`), so their filter chains cannot be read the
same way. Their filters are therefore unknown, not assumed equal.

## Reference Composition Audit: No Oral Decoy, And No Marker Under A chrY Fix Patch (2026-07-26)

Prompted by the observation that decoy/alt handling in a reference is designed around ~100 bp reads,
while ours average ~43 bp. Audited `mapping/index/hg38p14DH3630O.fa` (3.48 Gb) directly.

**No oral decoy is present**, despite the trailing `O` in the build name. Of 26,552 non-`chr*`
contigs, 26,355 are HLA allele sequences (126 Mb); there is no eHOMD content. `METHODS.md`'s
recommendation — omit oral decoys for this dataset, since they produced sparse spiky hits after
human-reference filtering — was in fact followed. Composition: 25 primary, 320 `_alt`, 2,349
unplaced/random/decoy, 1 EBV, 26,355 HLA.

**The real short-read risk in this build is the alt and fix contigs, not decoys.** `bwa aln` is not
ALT-aware (recorded in `mapping/map_se_adna.sh`), so the `.alt` companion is never consulted and
those contigs act as plain extra sequence competing for reads. Fix patches are the worst case: they
are *corrections* to the primary and therefore near-identical to it, so a read from a patched region
maps equally well to both and lands at MAPQ 0. At 43 bp there is less sequence to break the tie than
the reference designers assumed.

chrY carries three of them. `chrY_KZ208924v1_fix` (209,722 bp) aligns to chrY in **four** blocks —
not, as its 9.1–21.8 Mb min-to-max span suggests, contiguously:

    9,107,609  – ~9,112,000      (~4 kb)
    21,578,942 – 21,739,558      (160.6 kb)
    21,741,438 – 21,747,867      (6.4 kb)
    21,805,281 – 21,828,660      (23.4 kb)

**No marker in any of our three sets falls inside any block.** Nearest approach is `L166`/`L167` at
**15,077 bp** past the end of the last block; every other marker is ≥396 kb away. `Z6219`, `PF3239`
and `Z6494` sit in the large gap between blocks 1 and 2 and are 4.4–6.1 Mb clear.

This is a negative result and it was worth obtaining rather than assuming: had `L166` sat 15 kb
earlier it would have been under a near-identical duplicate of itself, and the observed `0% MQ0` at
that site would have been impossible. It also explains why the empirical `pct_mq0` figures are clean
— there is no chrY fix-patch sequence anywhere near the sites this project calls.

### Measured: what the alt/fix contigs actually take (2026-07-26)

`bwa aln` never reads the `.alt` file, so alt and fix contigs compete for reads with no ALT-aware
MAPQ correction. Whether that matters is an empirical question, and the answer is in our own BAMs.
Read density on the chrY contig family (`samtools idxstats`, reads per kb):

| contig | length | `I14677` | `I14678` | vs chrY primary |
|---|---|---|---|---|
| `chrY` (primary) | 57.2 Mb | 2.59 | 2.88 | — |
| `chrY_KI270740v1_random` | 37 kb | 69.8 | 75.5 | **~27x** |
| `chrY_KN196487v1_fix` | 101 kb | 15.3 | 16.4 | **~5.7x** |
| `chrY_KZ208924v1_fix` | 210 kb | 2.89 | 3.31 | 1.12–1.15x |
| `chrY_KZ208923v1_fix` | 48 kb | 1.63 | 1.47 | 0.55x |

The pull-away is real but lands where it does no harm. The 27x contig is unplaced Yq satellite.
The 5.7x contig aligns to chrY 56,821,509–56,887,903, which is **inside the already-declared no-go
region** `56.69–56.88 Mb Yq_het/PAR2`. And `chrY_KZ208924v1_fix` — the only one whose blocks come
near a marker (15 kb from `L166`) — takes reads at **background density**, i.e. it is carrying its
own share rather than stealing.

**A gap this exposes in our QC.** `ylib.mapq_audit()` reports `pct_mq0` and `n_mq_top` *at the site*.
That detects ambiguous reads which **stayed**. It is structurally blind to reads that were **pulled
away** to another contig, because those never appear in a chrY pileup at all. The table above is the
missing measurement and it was typed at a prompt, not committed — it belongs in a tool.

## Sensitivity Test: What Fraction Of Reads Can Even Come Back (2026-07-26)

`annotate/y_mappability.py`, results in `results/mappability/y_marker_mappability.tsv`. Every read of
length L overlapping each marker is cut from GRCh38 and mapped back with the pipeline's own
`bwa aln` options; `frac_recovered` is the fraction returning to the right locus at MAPQ>=25. Reads
are exact reference copies, so a failure is pure mappability with no sequencing error or damage
mixed in. 22 markers x 4 read lengths x 4 references.

**Mean `frac_recovered`:**

| target | 35 bp | 45 bp | 60 bp | 90 bp |
|---|---|---|---|---|
| `working` (hg38p14DH3630O) | 0.797 | 0.912 | 0.951 | 0.956 |
| `noalt` (GRCh38 no-alt masked) | 0.797 | 0.912 | 0.951 | 0.956 |
| `hs37d5` (GRCh37) | 0.797 | 0.912 | 0.951 | 0.956 |
| `chm13` (T2T) | 0.783 | 0.907 | 0.951 | 0.956 |

**1. The custom reference costs nothing and gains nothing here.** `working` and `noalt` are identical
to three decimals at every length. The alt, decoy and HLA content neither steals reads from our
markers nor protects them. This is the controlled version of the `idxstats` measurement above, and it
agrees with it.

**2. GRCh37 and GRCh38 are indistinguishable at these sites.** `hs37d5` matches GRCh38 exactly at all
four lengths. The concern that the regions under investigation might have changed materially between
builds is **answered in the negative for mappability** — which is the property that governs whether a
read is recoverable at all. It does not speak to sequence-level edits, but no marker changed its
recoverable fraction by so much as 0.001.

**3. Read length, not reference, is what costs power.** Mean recovery falls from 0.956 at 90 bp to
0.797 at 35 bp. A nominal 10x library at 35 bp reads is really interrogating these sites at ~8x.
Thirteen markers are at 1.000 in every build at 45 bp (`P15`, `M3308`, `PF3147`, `PF3177`, `L91`,
`P287`, `PF3239`, `Z6219`, `Z6287`, `Z6516`, `S19530`, `FGC5721`, `L167`), so the loss is
concentrated, not diffuse.

**4. `FGC5687` is uncallable and is in a marker set we use.** It is one of the three markers in
`markers/Z6494_exclusion.txt`:

    L=35  rec=0.000   35/35 reads MQ0   26 off-target -> chrX:9, chr6:3, chr17_alt:2
    L=45  rec=0.000   45/45 reads MQ0   34 off-target -> chrX:23, chr6:5, chr2:3
    L=60  rec=0.000   45/45 reads MQ0   23 off-target -> chrX:23
    L=90  rec=0.067   21    reads MQ0   11 off-target -> chrX:11

Identical in all four references, so it is an intrinsic property of the site, not a build artifact.
`FGC5687` sits in X-homologous sequence and **cannot be genotyped at any aDNA read length**. Every
previous "no coverage at `FGC5687`" in this project was mislabelled: it was never a coverage
shortfall, it is a site the aligner cannot place a short read on. The `Z6494` exclusion therefore
rests on `Z6215` (0.778 at 45 bp) and `Z6494` (1.000 at 60 bp+, 0.486 at 35 bp) — two markers, not
three.

**5. CHM13 carries a different base at `Z6488`.** Our catalogue has `anc=T der=C`; CHM13's chrY
reports **`A`** at the discovered position — neither allele. Every GRCh38/GRCh37 target reports `T`
correctly. CHM13v2's chrY is HG002's chromosome, a different individual and a different haplogroup
from GRCh38's, so this is the ref/alt inversion trap the check exists for. Whether it is a genuine
third allele or a one-base offset from a nearby indel has **not** been determined; either way,
**CHM13 must not be used for calling at `Z6488` until it is.** All other 21 markers match their
ancestral allele in all four builds.

**Method note.** The tool needs no liftover: the tile carries its own truth, so the modal implied
position is the marker's coordinate in the target build, discovered rather than asserted. This
yielded the GRCh38 -> T2T coordinates as a by-product (`L166` 21,843,737 -> 22,679,247; `PF3239`
15,205,748 -> 16,112,288; `Z6494` 17,131,187 -> 18,037,700) for an assembly we hold no chain for.

**One bug found and fixed before these numbers.** The first run reported `hs37d5` at 0.000 recovery
for every marker — hs37d5 writes `Y` where GRCh38 writes `chrY`, so every correctly-mapped read was
counted off-target. A tool bug indistinguishable from a catastrophic finding, which is the argument
for having a self-map control in the design: `working` against itself must return ~1.0, and any
target that returns 0.000 everywhere is a naming failure, not biology.

## The 1000G Phase 3 chrY Call Set Cannot Answer The Novel-Variant Question (2026-07-26)

The project's actual target is whether any public sample shares the Iceman's *novel* (uncatalogued)
derived variants, which is what would establish a shared sub-`L166` patriline. `HG02681`
(Pakistan/Punjab, YFull `G-Z6285`) is on this machine already, inside
`/mnt/GenomicData/1KG/1KGchrY/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz` (column 508
of 1,233 males), so the obvious question is whether it is worth testing.

It is not, and the reason generalises to every 1000G male.

Tested with `annotate/panel_membership.py --sites` (added for this) against the phase 3 chrY site
list (62,042 positions), lifting GRCh38 -> GRCh37:

- **Named markers: 17 of 22 present**, including `L166`, `L167`, `Z6219`, `Z6494`, `PF3239`. Absent:
  `FGC5696`, `FGC5721`, `FGC5687`, `Z6215`, `Z6488`.
- **Novel candidate positions: 0 of 21 present.** Not one.

So the call set can re-confirm known nodes and cannot address the question this project is actually
asking. Under the working data rule — public modern samples earn a place only if they yield genuine
insight — `HG02681` does not clear the bar, and no per-sample modern genotype needs to enter this
repository. That is a cleaner outcome than a judgement call about policy.

**Five of the 21 novel positions do not lift to GRCh37 at all:**

    chrY:10,768,171   10,964,462   10,990,649   10,996,925   11,667,647

All five fall in a ~1 Mb window at 10.7-11.7 Mb, and the last is inside the declared
`11.1-11.7 Mb` no-go region. "Does not lift" here means *no same-strand chain block*, which
`load_chain()` deliberately reports rather than guessing at: it may mean the sequence is absent from
GRCh37, or that it lies in an inverted block. Which one has not been determined.

**This is the one place the GRCh37/GRCh38 difference genuinely bites.** The mappability sweep found
GRCh37 and GRCh38 indistinguishable across all 22 *named* markers at every read length. That result
does not transfer here: a fifth of the novel candidate positions have no straightforward GRCh37
equivalent, so any analysis performed on GRCh37-mapped data — which is every depositor BAM in this
project, and `CGG017683` in particular — is structurally unable to see them. The sensitivity result
is scoped to catalogued markers in well-behaved sequence and must not be quoted as though it cleared
the novel scan.

## `Z6219` Read-Terminus Check: None Of It Is Damage (2026-07-26)

Tool: `annotate/y_read_evidence.py`. The registered rules score a single derived read at a `C>T`
site `nocall_damage_prone_1read` on **site class alone** — they never ask where in the read the
mismatch sits. Deamination is a terminal process, so that is the question that decides it.

**Every `Z6219`-derived read in the Oberbipp series, all six individuals:**

| sample | strand | read_len | BQ | MAPQ | dist_5p | dist_3p | NM |
|---|---|---|---|---|---|---|---|
| `MX210` | + | 68 | 41 | 37 | **30** | 37 | 1 |
| `MX210` | + | 31 | 41 | 37 | **28** | 2 | 1 |
| `MX210` | − | 44 | 41 | 37 | **15** | 28 | 1 |
| `MX210` | − | 50 | 41 | 37 | **34** | 15 | 1 |
| `MX187` | − | 53 | 41 | 37 | **22** | 30 | 1 |
| `MX211` | − | 76 | 41 | 37 | **72** | 3 | 1 |
| `MX212` | + | 39 | 41 | 37 | **20** | 18 | 1 |
| `MX213` | − | 41 | 41 | 37 | **14** | 26 | 1 |
| `SX10` | + | 51 | 41 | 37 | **46** | 4 | 1 |
| `SX10` | − | 59 | 41 | 37 | **48** | 10 | 1 |

`dist_5p` is distance along the **original molecule** — a reverse-strand read is stored
reverse-complemented, so its 5' end is at the right of the alignment, and naive
"distance from alignment start" would mislabel exactly the reads that matter.

Ten reads. The closest any of them comes to a 5' terminus is **14 bp**. Every one has `NM=1`, so
the site is the only mismatch on its molecule — no read shows the multi-mismatch pattern of a
genuinely damaged fragment. All are BQ 41 at the `bwa aln` MAPQ ceiling of 37.

**These libraries are UDG-treated, measured from the BAMs themselves:**

| sample | pos 0 | pos 1 | pos 2 | pos 3 |
|---|---|---|---|---|
| `MX210` | 0.0103 | 0.0008 | 0.0008 | 0.0005 |
| `MX213` | 0.0221 | 0.0014 | 0.0005 | 0.0006 |
| `SX10` | 0.0136 | 0.0011 | 0.0009 | 0.0006 |
| `MX187` | 0.0203 | 0.0011 | 0.0006 | 0.0001 |

C>T frequency by distance from the 5' end. Terminal damage survives at 1–2%; by position 1 it is
already an order of magnitude down and at background beyond. That is the partial-UDG signature, and
it is consistent across the study, as one lab's protocol should be. *(Rates are per read covering
the position, not per C — divide by the ~0.21 C content for a per-cytosine figure, so position 0 is
~5% per C and position 14+ is ~0.25%.)*

**Conclusion: the `Z6219` derived calls are real.** At 14–72 bp from the 5' end, in libraries whose
deamination is already at background by position 1, the chance that ten independent molecules across
six individuals all carry a damage-induced `C>T` is not worth computing. The four single-read calls
the registered rules refused were refused correctly on the information those rules use, and they are
now resolvable with information the rules do not look at.

**What this does to the topology.** Both constraints on `PF3239 → Z6219 → L166` are now solid:

- **Above `L166`** — Oberbipp: 10 derived reads at `Z6219` across four independent lineages
  (`MX210`, `MX213`, `SX10`, and Family A as one observation), none damage-positioned; ancestral
  11/11 at `L166`, a transversion and therefore damage-immune in both directions.
- **Below `PF3239`** — the Sardinians: `I14677` derived 5/5 at `PF3239` and ancestral 4/4 at
  `Z6219`; `I14678` 3/3 and 2/2. An *ancestral* call at a `C>T` site is damage-robust by direction,
  since deamination cannot manufacture the ancestral `C`.

**It is still post-hoc.** Nothing above was pre-registered, and the strength of the evidence is not
a substitute for having said in advance what would count. The pre-registration should now be written
knowing this, and should state what would falsify `PF3239 → Z6219 → L166` rather than what would
confirm it — the failure mode of the last one was assuming the labels it set out to test.

## Allele-Specific Reference Bias: An aDNA-Only Artifact Invisible To Modern Curators (2026-07-26)

Prompted by the observation that FTDNA omits `Z6219` and YFull rates it 3/5. That turned out to be
a false lead — YFull's own FAQ states that **2–5 stars are all "good quality" and all used to build
the YTree**, only 1 star being low — so 3 stars is unremarkable, and `Z6219` is clean under every
test below. The investigation it prompted found a real problem elsewhere.

**A limitation of `results/mappability/` is corrected here.** That sweep tiled reads cut from the
reference, so every read carried the **ancestral** allele. It measured ancestral-read recovery only.
`annotate/y_mappability.py --allele der` substitutes the derived base before tiling.

**Derived vs ancestral recovery, `noalt` reference, by read length:**

| marker | 45 bp | 100 bp | 150 bp |
|---|---|---|---|
| `Z6219` | +0.000 | +0.000 | +0.000 |
| `L166` | **−0.800** | +0.000 | +0.000 |
| `L167` | **−0.778** | +0.000 | +0.000 |
| `Z6208` | **−0.889** | −0.320 | +0.000 |

At 45 bp a derived read at `L166` is recovered 15.6% of the time against 95.6% for an ancestral
read. Identical in `working`, `noalt` and `hs37d5`, so it is not the custom reference.

**Mechanism — it is our own threshold, not paralogy.** Single-read diagnostic:

    L166_anc    MAPQ=37  XT=U X0=1 X1=0 NM=0   PASS
    L166_der    MAPQ=23  XT=U X0=1 X1=1 NM=1   FILTERED
    Z6208_der   MAPQ=20  XT=U X0=1 X1=2 NM=1   FILTERED
    Z6219_der   MAPQ=37  XT=U X0=1 X1=0 NM=1   PASS

Every alignment is **unique** — `XT:A:U`, `X0=1`. The derived reads are not ambiguous; they map to
one place, the correct place. MAPQ falls only because a *suboptimal* hit exists elsewhere (`X1>=1`)
and the read carries its single mismatch. `bwa aln`'s MAPQ conflates "uniquely placed" with "matches
perfectly", and our `MAPQ >= 25` cut therefore discards uniquely-mapped derived reads while keeping
every ancestral one. **The bias is allele-specific and its direction is toward calling ancestral.**

**Why no modern resource can warn about this.** The effect vanishes by 100 bp for `L166`/`L167` and
by 150 bp for `Z6208`. YFull's ratings are derived from 100–150 bp files and FTDNA's from Big Y; in
that regime these markers are flawless. This is the **second** independent way modern reference
resources fail to anticipate ancient data, alongside equivalence blocks that can only be split by a
sample on an intermediate branch (see the `Z6219` node work). Ancient Y calls are routinely checked
against trees built in a regime where the relevant failure modes do not occur.

**F6 checked and does not obtain — the Oberbipp result survives.** The obvious danger was that the
`L166` ancestral calls underpinning the refutation of the "7 of 10" claim were an artifact of
filtered derived reads. Registered in `PREREG_Z6219_node.md` Amendment 2 before being run, then
re-examined at `chrY:21843737` across all 15 Swiss BAMs **with no MAPQ floor**:

    14 reads, 8 samples (MX182, MX183, MX187, MX204, MX210, MX211, MX212, MX213)
    all base C (ancestral), all MAPQ 37, all NM=0, zero derived reads at any MAPQ

Nothing was being filtered. The ancestral calls stand as recorded and are now better established,
having survived a test designed to break them.

**What must change in interpretation, regardless.** These are arithmetic consequences, not findings:

- Every `L166` and `Z6208` **derived** count in this repository is a **lower bound**. `I5118`'s
  single derived `L166` read implies of the order of six underlying molecules — the evidence there
  is stronger than reported, not weaker.
- **`no_coverage` at `L166` or `Z6208` may not be read as evidence of anything.** A genuinely
  derived man is expected to show few or no reads at those sites. Several Oberbipp samples are
  `no_coverage` at `L166` and must not be counted toward the ancestral tally, which they never were.
- **Ancestral** calls are unaffected in principle: a derived molecule cannot produce an ancestral
  read, so the filter cannot manufacture the ancestral direction.

**The repair is identified and deliberately not applied.** Filtering on uniqueness (`XT:A:U` and
`X0=1`) rather than raw MAPQ is allele-neutral and is very likely the correct fix. Adopting it now,
having just discovered it changes specific results, is the post-hoc rescue this project's protocol
exists to prevent. It is recorded as a proposal requiring its own pre-registration and a re-run of
every affected dataset. Note also the tension named at the time: any such change — to the filter or
to the reference — recovers real signal at the cost of comparability with every other analysis.

## YFull's `(H)` Flag Independently Reproduced, And Test B's Real Power (2026-07-26)

YFull's live tree defines `G-L166` with **32 SNPs**; the list was supplied and is committed as
`markers/yfull_L166_defining.txt`. `Z6219` is **26th** in it, confirming the rank. Eleven of the 40
names (counting synonyms) do not resolve in our catalogue — all `Y######` YFull-assigned names plus
`Z6288` — because YFull publishes no coordinate or flanking sequence for SNPs it names itself. That
is a hard power limit on the registered Test B, not an omission.

**Allele-aware recovery across the 22 testable positions, `noalt`, 45 bp:**

| marker | pos (hg38) | anc | der | delta | note |
|---|---|---|---|---|---|
| `Z6220` | 15,115,629 | **0.000** | **0.000** | — | -> chr6, chr5, chr10 |
| `Z6519` | 21,049,898 | **0.000** | **0.000** | — | **YFull `(H)`** -> chr14:19, chrX:15 |
| `C101821` | 9,575,240 | **0.067** | **0.067** | — | -> chrX, chr4, chr5 |
| `Z6206` | 12,796,285 | 0.733 | 0.556 | −0.177 | |
| `Z6209` | 14,364,857 | 0.822 | 0.667 | −0.155 | |
| `Z6208` | 13,776,249 | 0.933 | **0.044** | **−0.889** | allele-biased |
| `L166` | 21,843,737 | 0.956 | **0.156** | **−0.800** | allele-biased |
| `L167` | 21,843,756 | 1.000 | **0.222** | **−0.778** | allele-biased |
| `Z1370` | 7,943,467 | 1.000 | **0.333** | **−0.667** | allele-biased |
| `Z6219` | 13,782,251 | 1.000 | 1.000 | +0.000 | clean |

Thirteen further markers (`FGC5671`, `FGC5672`, `FGC5675`, `FGC5721`, `FT91632`, `FT191098`,
`S10301`, `S19530`, `Z6134`, `Z6287`, `FGC5696`, `FGC5712`, `Z6499`) are at 0.867–1.000 in both
alleles.

**The `(H)` flag is independently reproduced.** `Z6519` is the only `(H)`-marked SNP resolvable in
our catalogue, and it comes back at **0.000 recovery in both alleles**, with reads scattering to
chr14, chrX and chr8. Our tool had no knowledge of the flag. That is external validation of
`y_mappability.py` against a curator's own annotation, and the first independent check this tool has
had.

**But the flag does not cover everything that matters here.** `Z6220` (0.000) and `C101821` (0.067)
are equally unusable and are **not** `(H)`-marked. `(H)` marks homology as it manifests in the
modern 100–150 bp regime; at 45 bp more sites fail than that flag anticipates. Same structural point
as the allele bias: curator annotations are produced in a regime where the aDNA failure modes do not
occur.

**`Z1370` joins the allele-biased set.** `L166`, `L167`, `Z6208` and now `Z1370` all lose most
derived reads at 45 bp. `Z1370` is not in any marker set this project currently uses; it is recorded
so it is not adopted later without the caveat.

**Test B's real power, fixed before the scan runs.** Of 32 SNPs defining the node: 11 untestable
(no public coordinates), 3 unusable at any allele (`Z6220`, `Z6519`, `C101821`), 4 severely
allele-biased, and **15 clean**. The co-segregation scan therefore has ~15 usable positions rather
than the ~11 estimated in the prereg, but the F3 threshold must be evaluated against 15, not 32 —
absence of co-segregation among markers that cannot be called is not evidence.

**Reproducibility check passed.** After making the filter parameters injectable, re-running
`y_genotype_batch.sh` at defaults reproduced **all five** `results/swiss15/` tables byte-identically,
and wrote `swiss_params.txt` recording thresholds, reference, inputs and git commit `44c6bdd` with a
clean tree.

---

## Tests A and B of `PREREG_Z6219_node.md`, run 2026-07-26

Both registered tests were run in the order A then B. Neither falsifier fired, and the co-segregation
scan returned one additional marker — but the power is thin enough that the qualifier matters more
than the headline.

### Test A — F1 does not obtain

`I5118` (Mezőcsát-Hörcsögös, Hungary, 3300–3000 BCE) is the **only `L166`-derived individual in
either cohort**, and he is derived at `Z6219` rather than ancestral:

| site | class | anc | der | registered call |
|---|---|---|---|---|
| `Z6219` chrY:13,782,251 | C>T | 0 | 1 | `nocall_damage_prone_1read` |
| `L166` chrY:21,843,737 | C>A | 0 | 1 | `DERIVED_1read_transversion` |

Both at MAPQ 37, 0% MQ0, `site_qc pass`, and on **different molecules** (`ERR2207344.6856395`,
`ERR2207344.6878821`). Across all 20 individuals in `results/testB/` and `results/testB_unhedged/`,
**no sample is `L166`-derived and `Z6219`-ancestral**. F1 had a real opportunity to fire and did not.

Per prereg §8 the C>T site requires read-terminus evidence, and it is decisive here: the derived read
sits **37 bp from the 5' terminus and 46 bp from the 3'** on an 84 bp molecule, while `I5118`'s
library runs 0.67% C>T at position 0 falling to **0.06% in the interior**. Damage probability for
that base is of the order of 1 in 1700.

**The registered call remains `nocall`.** One read is one molecule and the rule refuses it on site
class. The terminus evidence is reported beside the call, never folded into it.

**This is a survived falsification, not positive support.** H2 — `Z6219` genuinely `L166`-equivalent
— predicts the identical observation, because a man below both SNPs is derived at both under either
topology. Test A discriminates H1 from nothing; it only had the power to kill H1.

### Test B — F3 does not obtain, on one marker

Scan of the 22 testable positions of `markers/yfull_L166_defining.txt`, upstream cohort
`MX210, MX213, SX10, FamilyA` (pooled, one observation) against outgroup `I14677, I14678`:

| pattern | n | markers |
|---|---|---|
| `splits_with_Z6219` | 2 | `Z6219`, **`FGC5671`** |
| `stays_with_block` | 3 | `L166`, `L167`, `FGC5672` |
| `uninformative` | 17 | no coverage, single damage-prone read, or mixed |

**`FGC5671` (chrY:7,784,648, G>A) co-segregates with `Z6219`.** Upstream: `MX210` 2 derived,
`MX213` 1, `SX10` 1, `FamilyA` (`MX209`) 1 — **five derived reads, zero ancestral, four independent
lineages**. Outgroup: `I14678` 3 ancestral, `I14677` 1, and `I15942` 2 — six ancestral reads, zero
derived. This is P2 satisfied: the node does not rest on a single SNP.

Damage was checked and does not explain it, though the first look suggested it might. Four of the
five derived reads sit within 8 bp of the 3' terminus, which is where G>A deamination lives — but
**none is at position 0**, and these partial-UDG libraries collapse from 0.7–1.7% at position 0 to
0.1–0.2% at position 1 and beyond. Per-read damage probability is ≤0.2% for every one of the five.
Terminal *proximity* was the wrong statistic; terminal proximity weighted by the library's own rate
is the right one, which is exactly why prereg §8 requires both.

Mappability is clean at `FGC5671`: **1.000 recovery at both alleles** at 35, 45 and 60 bp, zero MQ0,
zero off-target — unlike `L166` (derived 0.000/0.156/0.867 at the same lengths). So neither the
upstream derived reads nor the outgroup ancestral calls are artifacts of allele-specific bias.

`MX210`'s two reads start 46 bp apart (7,784,587 and 7,784,633) despite consecutive read IDs, so
they are distinct molecules and F4 does not obtain for them.

### How weak this actually is

- **Only one registered `DERIVED` call.** Under the per-sample rules `MX210` alone is `DERIVED` at
  `FGC5671`; the other three lineages are `nocall_damage_prone_1read`. The co-segregation rests on
  one registered observation plus three consistent single-read nocalls.
- **All five upstream reads are + strand.** `y_markers_pileup.py` parses strand-agnostically and
  passes `1/1` to `site_qc`, so the `MARGINAL_single_strand` test cannot fire. Not adjudicated here;
  recorded because it is unflagged by construction.
- **17 of 22 positions are uninformative** — untested, not tested-and-negative.
- **F2 is largely untested, not passed.** The clause "ancestral at `PF3239`" could not be evaluated
  for `MX210`, `MX213` or `SX10`: all three have **no coverage** at `PF3239`. Only `FamilyA` is
  confirmed `PF3239`-derived (2 reads pooled, 0 ancestral). The "outside haplogroup G" clause does
  not obtain — all four lineages are derived at `P15`, `L91`, `Z6043`, `M3308`, `P287` with zero
  ancestral calls anywhere in the backbone.

### The gap that decides it, and which this scan cannot close

`FGC5671` derived in Oberbipp and ancestral in Sardinians is equally consistent with two things:
`FGC5671` sitting on the `Z6219` branch, **or** `FGC5671` being a private Oberbipp variant with no
relation to the `L166` block at all. Distinguishing them requires the marker in an **`L166`-derived**
man, and `I5118` — the only one available — has **no coverage at `FGC5671`**. Until that is closed,
"the node rests on two SNPs" is an inference from YFull's block assignment, not from our reads.

### Reproducibility

The Test B run regenerated every previously committed table as a side effect. **All ten** tables in
`results/swiss15/` and `results/unhedged/` reproduce byte-identically from
`results/testB/` and `results/testB_unhedged/`, at `git_commit=96de9d7`.

## chrY-only references manufacture confident false calls (2026-07-26)

The prediction recorded as UNTESTED in `RUNLOG.md` was run and holds.

| marker | whole-genome | chrY-only |
|---|---|---|
| **`FGC5687`** | **0.000** — 45/45 reads MQ0, 29 off-target | **1.000**, 0 MQ0, 0 off-target |
| `Z6215` | 0.778 | 1.000 |
| `Z6208` | 0.933 | 1.000 |
| nine others | unchanged | unchanged |

`FGC5687` is X-homologous. Mapped against a whole-genome reference its reads have somewhere else to
go, land ambiguously, and the site is correctly unusable. Mapped against chrY alone **every read
returns at full MAPQ with 0% MQ0**, because the paralogue has been removed from the reference.

**The failure is silent and self-certifying.** `pct_mq0` is the statistic that catches a collapsed
repeat, and here it reads 0% *because* the reference is deficient — the deficiency erases its own
evidence. A chrY-only pipeline emits a confident call with a clean QC flag at a position that cannot
be called at all. None of the three "improvements" above is real; each is the same artifact at a
different magnitude.

This is not only about our tooling. chrY-only mapping is common in ancient-DNA Y work, so published
calls at X-homologous positions can carry full confidence and no warning. It also bounds what our
own `Z6494` exclusion rests on: `FGC5687` returns 0.000 against every whole-genome reference tested,
so that exclusion rests on **two** markers, not three — and a Y-only reanalysis would appear to
recover the third while actually inventing it.

**Correction to an earlier note.** `RUNLOG.md` recorded that `bwa index` "did not complete under
load". It was killed when its parent shell exited, at 53% of the BWT after about two minutes. Under
`setsid` the same index builds in 199 s. Load was never the cause.

## E1: `FGC5671` is derived in the Iceman — and `Z6499` is not (2026-07-26)

Registered as E1 in `PREREG_Z6219_node.md` under `PROTOCOL_extending_analyses.md`, before any
Iceman read at either position was seen. Predicted outcome: derived at `FGC5671`. It is.

### The primary result

| | depth | ancestral | derived | pct_mq0 | site_qc | call |
|---|---|---|---|---|---|---|
| `FGC5671` chrY:7,784,648 G>A | 7 | 0 | 7 | 0% | pass | **DERIVED** |

Seven derived reads, four + strand and three −, lengths 59–137 bp, `NM=1` on every one with no
other C>T or G>A mismatch anywhere in the read. `FGC5671` is a G>A site, so deamination lives at
the **3'** end; the closest read sits 8 bp from it and the rest at 23, 33, 35, 48, 54 and 131 bp.
No read is at position 0 or 1. The Iceman library is UDG-treated and measures **0.10% C>T / 0.08%
G>A at position 0**, falling to 0.02% by position 2 — an order of magnitude below the Oberbipp
libraries. Seven independent molecules, both strands, at interior positions, in a library that
deaminates at 2 in 10,000 there: this is not damage.

Note that `y_read_evidence.py` labels each of these reads `derived(damage-pattern)`. That verdict is
assigned from the **site's mutation class**, not the read's position — it fires identically for the
read 131 bp from the 3' terminus. It flags what to check, not what was found.

**F3 does not obtain.** A second position now co-segregates with `Z6219` on read-level evidence in
three cohorts rather than two: derived in Oberbipp, ancestral in both Sardinians, and derived in a
man who is himself `L166`-derived. The disjunction recorded in "The gap that decides it" above is
resolved in favour of its first branch — `FGC5671` is a shared branch marker, not a private
Oberbipp variant. Homoplasy (H3) would now have to have struck twice, at two independent positions,
in the same lineages.

What this does **not** establish is that `FGC5671` sits beside `Z6219` *above* `L166`. The Iceman
is derived at `L166`, `Z6219` and `FGC5671` alike, so all three are on his path and he cannot order
them. He removes the "private variant" alternative; he does not place the marker.

### The secondary result, which was registered as the interesting one

Of the 22 testable positions in YFull's 32-SNP `L166` definition, the 9 already in
`markers/L166_defining.txt` are controls and returned **9/9 DERIVED at `site_qc pass`** — the run is
valid. Of the other 13, under the counting rule fixed in advance:

| class | n | markers |
|---|---|---|
| block-confirmed | 7 | `Z6134`, `FGC5671`, `Z6209`, `FGC5672`, `Z6135`, `Z6206`, `Z1370` |
| **block-refuted** | **1** | **`Z6499`** |
| untested | 5 | `Z6220`, `FT91632`, `FT191098`, `C101821`, `Z6519` |

All five untested positions fail the inherited 30% MQ0 rejection threshold (57–78% MQ0) or have no
coverage behind it. `Z6220`, `FT191098` and `Z6519` carry `DERIVED`-direction calls that the MAPQ
rule discards; they are **not** counted as confirmed, and the rule is not relaxed to admit them.

**`Z6499` (chrY:13,321,379, C>T): the Iceman is ancestral 10/0.** Five + strand and five −,
`NM=0` on every read, MAPQ 55–60, `site_qc pass`, 0% MQ0. An ancestral `C` at a C>T site cannot be
manufactured by deamination, so the call is damage-robust by direction.

Mappability was checked because a derived allele that cannot map produces exactly this picture.
`Z6499` is imperfect at short reads and asymmetric against the derived allele — 0.400 der / 0.429
anc at 35 bp, 0.778 / 0.867 at 45 bp, with off-target placements on chr3 — but **1.000 at both
alleles at 60 and 100 bp**. Six of the ten reads are 62–111 bp, in the stratum where recovery is
perfect and symmetric, and all six are ancestral. An 11% relative deficit at 45 bp cannot turn a
derived man into 0 derived reads out of 10. The call survives.

The other ancient data agree as far as they go: `MX213` is 1 ancestral and `I14677` 1 ancestral at
`Z6499`, which is why Test B scored it `uninformative`. Every ancient read anyone here has at this
position is ancestral. The Iceman is simply the first sample with the depth to make it a call.

**So YFull's `L166` definition contains at least one SNP that an `L166`-derived man does not carry.**
The economical reading is that `Z6499` sits **below** the `L166` node, on the branch leading to the
modern samples the block was built from, and the Iceman's lineage diverged above it.

That is the same phenomenon as the `Z6219` hypothesis, pointing the other way. `Z6219` would split
the block upward, `Z6499` splits it downward, and both are invisible to a modern-only tree for the
same structural reason: **an equivalence block can only be split by a sample on the intermediate
branch**, and ancient samples are where those branches survive. The `L166` block is being opened at
both ends by the same defect in the same direction of inference. This strengthens the general claim
in §3 of the pre-registration — that this is a systematic confounder for ancient `G-L166` calls —
while remaining silent on `Z6219`'s own placement, which `Z6499` does not bear on.

Registered in advance as the notable outcome, so it is not a post-hoc finding. What it is *not* is
a test of anything: one refuted marker out of 13 in one individual, and the direction of the
inference (below rather than parallel) rests on parsimony, not on a sample that carries `Z6499`
derived and `L166` ancestral. No such sample has been looked for.

## E2: `CGG017683` cannot be placed relative to `Z6219` — H0, as registered (2026-07-26)

Registered as E2 in `PREREG_Z6219_node.md` with the expected result stated in advance as H0 at
~70–97% depending on the outcome asked for. That is what it returned.

### The result

3 of 20 testable positions covered, all single-read — against 17.5% expected at 0.192x, so **15.0%
observed**, exactly on model.

| marker | b37 pos | class | reads | call |
|---|---|---|---|---|
| `FGC5721` | Y:18,392,335 | T>G transversion | 1 derived | **`DERIVED_1read_transversion`** (a registered call) |
| `Z6134` | Y:6,702,979 | C>T | 1 derived | `nocall_damage_prone_1read` |
| `S19530` | Y:17,320,353 | T>C | 1 derived | `nocall_1read_transition` |

**`Z6219`, `L166`, `L167` and `Z6499` are all uncovered.** Not thinly covered — **zero reads**, and
that was checked with every filter switched off (`samtools view -c`, and `samtools depth -q 0 -Q 0`,
both returning 0 at all four). This matters because the deposit ran `samtools calmd -Erb` without
`-A`, so a BAQ-capped base quality could have produced empty pileups that merely looked like absent
molecules. It did not: there are no molecules there to filter.

**Verdict: H0.** The registered no-power criterion — fewer than 2 reads at `Z6219` *and* no
registered transversion call at `L166`/`L167` — is met on both limbs. The one registered call,
`FGC5721` derived, is uninformative for the question: the Iceman is derived at `FGC5721` too, so it
sits on the shared path and orders nothing.

**This sample cannot answer whether it is on the intermediate branch, and no additional work on
public data will change that.** The limit is 54,685 chrY reads. The user's observation that YFull
hangs it at `G-L166` but outside the resolved sub-branches stands unexplained and untestable here;
"outside the sub-branches" at this depth is equally well explained by missing coverage at the
defining SNPs, which is the same demotion-toward-the-parent artifact already noted for `ERS257168`.

### Two things worth keeping

**The `S19530` derived read reproduces independently.** The 2026-07-25 screen reported this sample
derived at `S19530/Z6213` off a single read, via an ad-hoc hg38→hg19 route whose output was never
committed. It is recovered here by a different chain — `annotate/y_lift_markers.py` against the
UCSC chain file, then `y_sites_pileup.py` — landing on the same base. That is a positive control on
the new lift path, obtained for free.

**Two of the 22 positions do not exist in GRCh37.** `FT91632` (hg38 10,756,989) and `FT191098`
(hg38 10,874,130) fall in no same-strand chain block and are reported `unmapped` rather than given a
coordinate. This is the concrete form of the loss described under "Every Sample In This Project Came
From A Mapped-Only Deposit": sequence that exists only in GRCh38 cannot be represented in a GRCh37
deposit at all, so those positions are permanently untestable in this sample no matter how deep it
were sequenced. Both were also `untested` in the Iceman on the MQ0 rule, for unrelated reasons.

The other 20 lifted with the hs37d5 reference base matching the expected ancestral allele at
**every one**, which is the check that would have caught a mis-lift.

## E3: the uniqueness filter is rejected — it trades an ancestral bias for a paralogue (2026-07-26)

Registered in `PREREG_uniqueness_filter.md` with the adoption criterion fixed before the measurement.
**U1 fails and U2 fails. `MAPQ >= 25` stands.**

Measured on simulated tiles cut from the reference — exact copies, no error, no damage — so any
ancestral/derived difference is produced by the filter and nothing else. 22 markers × 4 read lengths
× 2 references × 2 alleles.

### U1 — the bias is reduced but not repaired

Of 176 marker/length/target cells, **28** are asymmetric beyond 0.05 under `MAPQ >= 25`. Under the
uniqueness criterion **18** still are. No marker's derived recovery got worse, so the direction of
the repair is right; the magnitude is not.

| marker | len | MAPQ anc/der | uniq anc/der |
|---|---|---|---|
| `L166` | 45 | 0.956 / **0.156** | 1.000 / **0.689** |
| `L167` | 45 | 1.000 / 0.222 | 1.000 / 0.689 |
| `Z6208` | 45 | 0.933 / **0.044** | 1.000 / 0.933 |
| `Z6208` | 100 | 1.000 / 0.680 | 1.000 / 1.000 |

`L166` at 45 bp reproduces Amendment 2's quoted 0.956 / 0.156 exactly, which is a reproducibility
check on the whole measurement obtained for free. Uniqueness lifts derived recovery from 0.156 to
0.689 — a real improvement, and still **six times** the registered 0.05 threshold. U1 fails.

### U2 — and it rescues reads that belong somewhere else

`FGC5687` was named in the registration as the control for exactly this, before the run:

| marker | len | MAPQ | uniq | off-target reads |
|---|---|---|---|---|
| `FGC5687` | 35 | 0.000 | 0.000 | 26 / 27 |
| `FGC5687` | 45 | 0.000 | 0.000 | 34 / 30 |
| **`FGC5687`** | **60** | **0.000** | **0.250** | 23 / 23 |
| **`FGC5687`** | **100** | **0.160** | **0.890** | 6 / 5 |

At 100 bp the uniqueness criterion turns a marker with 0.160 recovery into one with 0.890, at a
position where the whole-genome measurement puts a quarter of the tile on another contig. Ancestral
and derived move together (0.890 / 0.890), so this is not an allele effect at all — it is a wholesale
admission of reads whose alternative locus is real. 29 cells fail U2. **Hard fail, no override**, as
registered.

### The trade is exact, and that is what makes it fatal

Uniqueness repairs the `L166`/`L167`/`Z6208`/`Z6494` asymmetry **only at 60 bp and above**, where all
four go to 1.000 / 1.000. `FGC5687` admits the paralogue **only at 60 bp and above**. The read
lengths at which the criterion fixes the bias are precisely the read lengths at which it stops
distinguishing the locus from its homologue. There is no window where it is safe and useful.

This is the same defect as the chrY-only reference finding, arriving by a different route. There, the
paralogue was deleted from the *reference*, so reads had nowhere else to go and `pct_mq0` read 0%
because the reference was deficient. Here the paralogue is deleted from the *acceptance rule*: `X0`
counts best hits and ignores that a suboptimal hit exists, so `XT:A:U` reads clean at a position that
cannot be called. Both make a read look unique by removing the thing it competes with, and both
produce a confident call with clean QC. `FGC5687` is the marker that exposed the first defect and it
is the marker that exposes this one.

### Consequences

**`MAPQ >= 25` remains in force unchanged.** `ylib.uniqueness_audit()` continues to report without
acting. No call in this repository moves, and the tooling change is additive: `y_mappability.py`
gained `n_unique` and `frac_recovered_unique` as reported columns, and every existing column
computes exactly as before.

**Amendment 2's consequence is now permanent, not provisional.** Every `L166`/`Z6208` derived count
here is a **lower bound**, and no `no_coverage` at those positions is evidence of anything — and the
obvious repair has now been tested and rejected, so that caveat cannot be retired by adopting it.

**U3 was not run.** The registered rule is adoption iff all three hold; U1 and U2 both failed, so the
decision was already determined and running the known-answer regression under a criterion that will
not be used would have proved nothing.

**The effect on this project's own calls was deliberately not computed.** The registration promised
it would not enter the decision. With the criterion rejected, computing it would produce nothing but
a table of what could be gained by adopting a filter that fails its safety control — which is the
temptation the firewall exists to block. Not computing it is the safer choice and is recorded as a
choice, not an omission.

**`Z6219` is untouched either way.** It recovers 1.000 / 1.000 at 35, 45, 60 and 100 bp under both
criteria, with zero off-target placements. The central evidence of `PREREG_Z6219_node.md` does not
depend on this filter question at all — which is worth stating precisely because the outcome was not
known when the measurement was registered.

**What would actually be right** is a criterion that examines the *identity* of the suboptimal hit —
whether the alternative locus is a real paralogue or a chance similarity — which `X0`/`X1` counts
discard by construction. That is a different aligner's problem, not a threshold change, and nothing
here authorises attempting it.

## Registered outcome of `PREREG_Z6219_node.md`: H1 supported (2026-07-26)

The verdict §8 demanded had never been written down. It is now recorded in the pre-registration
itself: **H1 supported**, F1/F3/F4 do not obtain, **F2 untested**.

The substance in one paragraph: the ordering `PF3239 → Z6219 → L166` rested on a single SNP, which
is precisely what recurrence looks like. It now rests on two. `FGC5671` is ancestral in the
Sardinians (below `PF3239`) and derived in Oberbipp men who are `L166`-ancestral (above `L166`), so
it sits on the same intermediate branch as `Z6219`; E1 ruled out the alternative that it was a
private Oberbipp variant by finding it derived in the Iceman. Homoplasy would have to have struck
twice at two independent positions in the same lineages.

**The registered vocabulary is coarser than the evidence and was not adjusted to fit.** All four
lineages are from one site, `FGC5671`'s upstream support is one registered call plus three
single-read nocalls, and F2 could not be evaluated for three of the four. "H1 supported" names which
of four boxes this falls in; it is not a claim that the node is established. The test that would
discriminate H1 from H3 is a second population, registered as P3 and not run.

## Test C is Registered, and the Power Analysis Changed What It Is (2026-07-26)

`PREREG_testC_aesch_muttenz.md`, written before any Aesch or Muttenz read was staged. Tables in
`results/testC_power/`. Test C is the only registered test that separates H1 (a real node between
`PF3239` and `L166`) from H3 (homoplasy, or a variant private to the Oberbipp patriline), because
everything supporting the current `H1 supported` verdict is one patriline at one site.

**The cohort named in the earlier documents does not exist.** "The 21 hedged Aesch/Muttenz
individuals", used in `PREREG_Z6219_node.md` §7 and in `RUNLOG.md`, conflates two sets. There are
exactly 21 individuals labelled `G-L166*` in the compilation, but **6 of them are the Oberbipp men
already genotyped in `results/swiss15/`**. The untested cohort is **13 Aesch + 2 Muttenz = 15**.

**Aesch has Oberbipp's kinship problem, and that is the real result of the power analysis.** The 13
Aesch candidates are not 13 observations: six of them (`Aes1`, `Aes12`, `Aes19`, `Aes20`, `Aes21`,
`Aes23`) are one documented family, and the rest fall into three more. What Test C buys is **seven
kin groups across two sites**, against the current verdict's one site. The count of individuals was
never the interesting number.

**Only 5 of the 22 YFull `L166`-defining positions are testable in this data at all.** Across all 15
Oberbipp libraries — 135,398 chrY reads at MAPQ >= 25 — `L166` attracted 14 reads, `Z6219` 10,
`L167` 7, `FGC5671` 5, `FGC5672` 3, and **ten positions attracted zero**. Four of those ten
(`FGC5696`, `FGC5721`, `Z6516`, `S19530`) *are* on the 1240k panel. So panel membership does not
predict coverage here, and the idea of buying power by pooling more markers is not available. The
testable set is `L166`, `L167`, `Z6219`.

**Registered expectation, for later comparison against the outcome:** ~4.0 individuals callable at
`Z6219`, ~5.9 at `L166`, and **~2.1 callable at both** — the last being the number that matters,
since the intermediate-branch signature is derived at `Z6219` and ancestral at `L166` *in the same
man*. Power to observe the split P3 predicts is 0.79 at a derived frequency of 0.5 and 0.33 at 0.1
or 0.9.

**The dominant uncertainty is the marker rate, and it is large.** The `Z6219` rate is fitted on ten
reads; its Poisson interval spans a factor of 3.8, which puts expected callable between **1.3 and
7.9** and split power between **0.23 and 0.98**. Test C is somewhere between barely worth running
and decisive, and which one cannot be known until the reads are mapped. It is registered anyway,
because the alternative is leaving the verdict resting on a single site.

**The model reproduces the cohort it was fitted on.** Run against Oberbipp/Rapperswil themselves it
predicts 2.08 callable at `Z6219` (observed 2) and 3.15 at `L166` (observed 4).

**One estimator choice worth recording, because it was nearly a silent bias.** The coverage proxy is
the compilation's `NRY` column, and `MX182`'s row carries `NRY = 2` against 9,148 mapped chrY reads
— a defect already noted here on 2026-07-25. Under a through-origin fit that single row raised the
proxy ratio from 0.6514 to 0.7455, inflating every power number by 14%, because a sample with a
near-zero denominator and real reads dominates the numerator. `y_power_estimate.py` now uses the
median of per-sample ratios and reports the through-origin value beside it. `MX182` is not dropped;
it stays visible in the committed table and simply stops dominating the fit.

**What this test may not claim, registered in advance.** The cohort is defined by the compilation's
`G-L166*` label, which is a hobbyist re-derivation, **not** by the publication's `PF3239` label that
P3 named. Furtwängler 2020 Supplementary Table 1 has not been re-read for Aesch, so the overlap
between the two sets is unknown and is listed as a required pre-staging check. Confirming or
refuting the compilation's labels is not a finding about the haplotree; the labels are a sampling
frame, chosen because it is reproducible from a committed file.

### Test C §7 pre-staging checks, against the publication itself (2026-07-26)

Furtwängler 2020's corrected Supplementary Information fetched and read directly (25 pp.,
Supplementary Tables 1 and 5) rather than via the compilation or a search summary. **The cohort
survives all three checks and nothing in the registration is amended.** Four things came out of it.

**1. The ID-mapping hazard does not exist for this cohort, and the assumption was backwards.**
Supplementary Table 1 uses `Aesch1`, `Aesch12`, `RA61` — *identical to the ENA `sample_alias`
values*. It is the compilation that abbreviates to `Aes1`, `Aes12`. So the `Aes`→`Aesch` mapping
recovers the publication's own IDs and the ENA↔publication correspondence needs no inference.

**2. `RA61` and `RA62` are male.** Supplementary Table 1 is titled "Y chromosomal haplogroup
assignment for **all male individuals**" and both are in it. The empty `Sex` field is a gap in the
compilation, not an ambiguity in the source.

**3. The registered cohort is a strict superset of the set P3 named.** Twelve of the 15 candidates
are published `PF3239` — which is *all ten* at Aesch and *both* at Muttenz, confirming the per-site
counts recorded above. The other three are published shallower: `Aesch6` `PF3147`, `Aesch7` and
`Aesch20` `FGC7739/Z6488`. They are retained, because they were selected by a stated reproducible
criterion and dropping them on a label read afterwards is rule-shopping; `Aesch6`'s predicted power
is 0.036, so retaining it costs nothing. Test C can now report under either cohort definition.

**4. The "seven kin groups" figure is the compilation's, not the publication's — and this is the
one that needed catching.** Supplementary Table 5 gives exactly one male–male first-degree pair
among the candidates, `Aes12`–`Aes19` (`pi_HAT` 0.427, 272,655 SNPs, `same Y HG = yes`), which is
precisely and only the merge the registration had already fixed. It contains **no entry pairing
`Aes1`, `Aes20`, `Aes21` or `Aes23` with anything**. The compilation's twelve-member "Family D" is a
re-derivation of unstated provenance, like its `Y-Haplotree-Variant` calls.

The correction runs in the cohort's *favour* — on the publication's evidence there are 14 units with
no documented relationship rather than seven kin groups — and it is deliberately not banked. Table 5
is scoped to first and second degree, so absence from it is not evidence of unrelatedness, and
neither source measures patrilineal sharing beyond third degree, which is the thing that actually
matters. FC3 stands.

**Two observations recorded in passing.** `Aes12`–`Aes24` are first-degree with `same Y HG = same
clade`, yet published `PF3239` and `Z6488` respectively — two men who must carry the same Y
chromosome, labelled one node apart. That is the "published terminal is a coverage floor, not a
placement" point of Table 5, now seen at Aesch as well as Oberbipp. Separately, **`SX11` is
internally inconsistent in Table 1**: terminal mutation `PF3239` against `YHG = G2a2b2a1a1`, and
`PF3239` is `G2a2a1a2a1`. `SX11` is the single Niederried individual in the published `PF3239` tally
of 20, so that tally should treat it as unresolved.

### The coverage proxy failed a second time, and this time in the cohort's favour (2026-07-26)

Recorded while mapping was still running, **before any Test C genotype existed**, because a coverage
observation noted before the alleles are read means something a coverage observation noted afterwards
does not.

The first two Test C libraries to finish give the first out-of-sample test of the proxy calibration
in `PREREG_testC_aesch_muttenz.md` §4. Predicted chrY reads are the compilation's `NRY` value times
the fitted ratio 0.6514; observed are reads on `chrY` at `MAPQ >= 25` in the `.rmdup.bam`:

| sample | `NRY` proxy | predicted chrY | observed chrY @ MQ25 | obs/pred |
|---|---|---|---|---|
| `Aesch12` | 20,241 | 13,185 | 12,082 | 0.92 |
| `Aesch13` | 17 | 11 | 15,953 | **1441** |

`Aesch12` validates the model — 8% under, well inside the Poisson interval that dominates every
number in §4.

`Aesch13` does not, and the failure is the *same defect* as `MX182`: a near-zero `NRY` against a
library with real chrY data. `MX182` carries `NRY = 2` against 9,148 mapped chrY reads; `Aesch13`
carries `NRY = 17` against 15,953 at `MAPQ >= 25`. The MAPQ filter removes only 3% of its chrY reads,
so this is not multi-mapping junk inflating a count — it is a normal library that the proxy column
reports as empty. Two instances in the only two cohorts where the column has been checked against
mapped reads is not a stray typo; it is a column that cannot be trusted per-row.

**§4 is not rewritten.** The power statement was fixed before staging and its numbers stand as what
was predicted. What changes is only the interpretation: `expected_callable = 4.0 of 15` at `Z6219` is
now known to be an *under*-estimate, because one candidate registered at essentially zero power has a
chrY yield in the cohort's upper half. Re-fitting the proxy now and quoting the larger number as the
registered power would be exactly the rule change after seeing its effect that
`PROTOCOL_extending_analyses.md` prohibits.

`PREREG_testC_aesch_muttenz.md` §7 check 4 required `Aesch13` to be reported and not quietly dropped,
on the expectation that it would return no power. That rule was written to stop a null being hidden.
It has instead retained one of the cohort's better libraries — which is the argument for the rule,
not against it: the same discipline that forbids dropping an inconvenient sample forbids dropping a
convenient one.

Whether those reads land on the decisive positions is a separate question, and the genotyping answers
it. chrY yield is necessary for power, not sufficient — ten of the 22 YFull `L166`-defining positions
attracted zero reads across all 15 Oberbipp libraries, four of them despite being on the 1240k panel.
