# Run Log (Otzi Reanalysis)

This file is the operator-facing ledger for reproducibility. Keep methods in `METHODS.md`, interpretation in `NOTES.md`, and put exact run invocations/results here.

## Scope
- Dataset: Otzi 2023 resequencing BAMs, EAGER-style prefixed reads (`M_/F_/R_`), published in filtered form.
- Goal: reconstruct PE semantics where possible, remap, evaluate duplicate handling, and compare variant-calling branches.

## Privacy Guardrail
- This file is for ancient/public runs only.
- Do not record modern/private sample IDs, direct file paths, or sample-specific call summaries here.
- For modern/private runs, store outputs outside the repo and keep only generalized commands/templates in tracked docs.

## Canonical Inputs
- Baseline BAM: `<raw_data_dir>/iceman.oetzi.UDG_D2049_combined.mapped_rmdup.bam`
- Reference (primary run): `mapping/index/hg38p14DH3630O.fa`

## Path Conventions
- `<raw_data_dir>`: original input BAM location.
- `<analysis_dir>`: run-specific BAM outputs used as DeepVariant input.
- `<legacy_data_dir>`: legacy SOLiD-era VCF/marker files.
- `<output_dir>`: final DeepVariant outputs and comparison reports.
- `<repo_root>`: checked-out path to this repository.

## Key Decisions Locked
- Keep `mapping/revert-bam.sh` stage-gated (`ENABLE_ALIGN`, `ENABLE_MARKDUP`).
- Use explicit `samtools fastq` split stream before `bwa mem -p` for mixed PE/SE handling.
- For DeDup, use primary-assembly slice for production branch.
- Treat DeDup output as one branch; singleton-heavy over-pruning is a documented caveat.
- For Y dual-liftover experiments, default to non-extended refs; only opt into H3630/decoy references intentionally.
- Treat compressed-reference dictionary naming (`ref.dict` vs `ref.fa.dict`) as a required preflight item.

## Confirmed Counts
### D2049 prefix counts
- Source BAM: `M=138980507, F=5545558, R=3582337`
- `pair.prim.bam`: `M=136687692, F=5138191, R=3238303`
- `pair.prim_rmdup.bam`: `M=131715629, F=3743710, R=2869582`

### D2049 DeDup removals (log-consistent)
- Merged removed: `4972063`
- Forward removed: `1394481`
- Reverse removed: `368721`
- Total removed: `6735265`
- Exact rate: `4.64%` (`6735265 / 145035948`)

### Final merged pass (D2049-D2052)
- Before final DeDup: `668715804`
- After final DeDup: `656719414`
- Removed final pass: `11996390` (`1.79%`)

## DeepVariant Status
- Ancient DNA branch: long `postprocess_variants` tail observed.
- Evidence of progress: active `pread64` streaming + high CPU in protobuf `_message.so`.
- Runtime caveat: non-primary contig scope with gVCF can create heavy postprocess overhead even with zero-read contigs.

## Next Commands (shortlist)
1. Finalize Ancient DNA DV outputs, then snapshot stats:
```bash
bcftools stats <output_dir>/iceman.vcf > <output_dir>/iceman.vcf.stats
```
2. Run comparison branch without final merged DeDup (same DV settings), then diff callsets.
3. WES/WGS modern runs: enforce `--regions` on primary assembly; keep shard count conservative on laptop.

## DeepVariant A/B Plan (small_model on/off)
Run both branches on the same BAM/reference and same region list.

Region syntax note:
- `run_deepvariant` accepts space-separated contig list after `--regions` (confirmed in current runs).

Primary contig list used:
```bash
REGIONS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
```

> **Correction (2026-07-26): the two commands below were never executed as written.**
> Reconstructed from `~/.bash_history` and the DeepVariant log directories, the runs that
> actually produced the files on disk were:
>
> | branch | `--output_vcf` | `--disable_small_model` | `--regions` | logs |
> |---|---|---|---|---|
> | A | `/output/iceman.vcf` | `false` | *(none)* | `~/logs/` |
> | B | `/output/iceman-nosmall.vcf.gz` | `true` | `"chr1 ..."` | `~/logs-nosmall/` |
>
> `iceman_sm_on.vcf.gz` / `iceman_sm_off.vcf.gz` never existed. Both runs used
> `...sort_rmdup.coord.bam`. DeepVariant wrote BGZF into the bare `.vcf` name for branch A, which
> is why `~/iceman.vcf` is compressed and unindexed while `~/iceman-nosmall.vcf.gz` has a `.tbi`.
> A third invocation targeted `iceman.vcf` with `--disable_small_model=true`; it aborted
> (`~/logs2/call_variants.log` is 0 bytes, timestamped ~2 h *after* `iceman.vcf` was written) and
> did not overwrite it. Identity of the surviving files is confirmed independently by chrY hom-alt
> `GQ` max (58 = small model on, 25 = off) and by the `438` count at `PASS,GT=1/1,DP>=5,GQ>=30`.
>
> The lesson is documentation drift, not renaming: what was committed here was the *intended*
> command, and the executed one differed. Record commands by capturing them at run time.

Branch A (`small_model` enabled):
```bash
/snap/docker/3377/bin/docker run --rm --user "$(id -u):$(id -g)" \
  -v <analysis_dir>:/input \
  -v <repo_root>:/repo \
  -v <output_dir>:/output \
  google/deepvariant:1.10.0 /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/repo/mapping/index/hg38p14DH3630O.fa \
  --reads=/input/iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam \
  --output_vcf=/output/iceman_sm_on.vcf.gz \
  --output_gvcf=/output/iceman_sm_on.g.vcf.gz \
  --num_shards=8 \
  --vcf_stats_report=true \
  --disable_small_model=false \
  --logging_dir=/output/logs_sm_on \
  --haploid_contigs=chrX,chrY \
  --par_regions_bed=/input/GRCh38_PAR.bed \
  --intermediate_results_dir=/output/dv_tmp_sm_on \
  --regions $REGIONS
```

Branch B (`small_model` disabled):
```bash
/snap/docker/3377/bin/docker run --rm --user "$(id -u):$(id -g)" \
  -v <analysis_dir>:/input \
  -v <repo_root>:/repo \
  -v <output_dir>:/output \
  google/deepvariant:1.10.0 /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/repo/mapping/index/hg38p14DH3630O.fa \
  --reads=/input/iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam \
  --output_vcf=/output/iceman_sm_off.vcf.gz \
  --output_gvcf=/output/iceman_sm_off.g.vcf.gz \
  --num_shards=8 \
  --vcf_stats_report=true \
  --disable_small_model=true \
  --logging_dir=/output/logs_sm_off \
  --haploid_contigs=chrX,chrY \
  --par_regions_bed=/input/GRCh38_PAR.bed \
  --intermediate_results_dir=/output/dv_tmp_sm_off \
  --regions $REGIONS
```

Minimal comparison commands:
```bash
bcftools stats <output_dir>/iceman_sm_on.vcf.gz  > <output_dir>/iceman_sm_on.vcf.stats
bcftools stats <output_dir>/iceman_sm_off.vcf.gz > <output_dir>/iceman_sm_off.vcf.stats

bcftools isec -p <output_dir>/iceman_sm_cmp \
  <output_dir>/iceman_sm_on.vcf.gz <output_dir>/iceman_sm_off.vcf.gz

wc -l <output_dir>/iceman_sm_cmp/0000.vcf \
      <output_dir>/iceman_sm_cmp/0001.vcf \
      <output_dir>/iceman_sm_cmp/0002.vcf \
      <output_dir>/iceman_sm_cmp/0003.vcf
```

### Implemented VCF comparator (small_model on/off)
Reusable script:
- `annotate/compare_vcf_runs.sh`

Example command:
```bash
cat > /tmp/iceman_y_markers.tsv <<'EOF'
chrY	12915617	M201
chrY	13776249	Z6208
chrY	19483669	L91
chrY	21843737	L166
EOF

annotate/compare_vcf_runs.sh \
  -a <output_dir>/iceman.vcf \
  -b <output_dir>/iceman-nosmall.vcf.gz \
  --label-a small \
  --label-b nosmall \
  --out-dir /tmp/iceman_smallmodel_compare \
  --regions-file /tmp/iceman_y_markers.tsv
```

Observed Iceman A/B results:
- all sites:
  - totals: `5011174 / 5011174`
  - shared: `5010771`
  - private: `403 / 403`
  - shared GT diff: `38024` (`0.758845%`)
- PASS-only:
  - totals: `3987840 / 3970147`
  - shared: `3969848`
  - private: `17992 / 299`
  - shared GT diff: `671` (`0.016902%`)
- Y marker check (`M201`, `Z6208`, `L91`, `L166`):
  - all four are `PASS` and `GT=1/1` in both runs.
  - notable confidence shift at `Z6208`: `GQ 39 -> 19`, `QUAL 38.5 -> 18.9`.

Saved result artifacts:
- `results/iceman_smallmodel_compare.summary.tsv`
- `results/iceman_smallmodel_compare.y_markers.tsv`

## Y Haplogroup Comparison Run
### Input branches compared
- Legacy SOLiD-era test VCF:
  - `<legacy_data_dir>/chrY_raw_Iceman_tst_hg38.vcf.gz`
- Illumina reanalysis DeepVariant VCF:
  - `<output_dir>/iceman.vcf` (gzip-compressed VCF without `.gz` suffix)
- Marker panel:
  - `<legacy_data_dir>/snps_hg38.vcf.gz`

### Commands used
```bash
# SOLiD branch (strict-quality pass used in comparison notes)
annotate/y_haplo_from_vcf.sh \
  -i <legacy_data_dir>/chrY_raw_Iceman_tst_hg38.vcf.gz \
  -o /tmp/iceman_solid_strict \
  --ann-vcf <legacy_data_dir>/snps_hg38.vcf.gz \
  --site-filter-mode any \
  --min-gq 30 --min-dp 10

# Illumina DeepVariant branch
annotate/y_haplo_from_vcf.sh \
  -i <output_dir>/iceman.vcf \
  -o /tmp/iceman_illum \
  --ann-vcf <legacy_data_dir>/snps_hg38.vcf.gz \
  --site-filter-mode deepvariant \
  --min-gq 20 --min-dp 3

# Candidate clade scoring across branches
annotate/y_clade_consistency.py \
  --input solid=/tmp/iceman_solid_strict.derived.tsv \
  --input illum=/tmp/iceman_illum.derived.tsv \
  --candidate G \
  --candidate G2a \
  --candidate G2a2a1a2a1a1b \
  --candidate G-Z6208 \
  --out /tmp/iceman_clade_compare.tsv
```

### Key outcomes recorded
- Illumina branch key marker states:
  - `M201` derived (`GT=1/1`, `GQ=42`, `DP=10`)
  - `L91` derived (`GT=1/1`, `GQ=23`, `DP=7`)
  - `L166` derived (`GT=1/1`, `GQ=14`, `DP=3`)
  - `Z6208` derived (`GT=1/1`, `GQ=39`, `DP=9`)
- SOLiD branch:
  - corresponding markers unresolved/no-call in current branch outputs.
- `y_clade_consistency.py` (from `/tmp/iceman_clade_compare.tsv`):
  - `illum` shows positive net support for `G`, `G2a`, `G2a2a1a2a1a1b`, `G-Z6208`.
  - `solid` shows strongly negative net support for same candidates under current filter regime.

### Marker-state/Path scoring extension (VCF+gVCF aware)
Added Python tools to keep both VCF and gVCF inputs supported while preserving `ancestral` and `nocall` evidence:
- `annotate/y_haplo_from_markers.py`
  - outputs `<prefix>.marker_status.tsv` with per-marker status:
    - `derived|ancestral|nocall|ambiguous`
  - outputs `<prefix>.derived.tsv` and `<prefix>.summary.txt`
- `annotate/y_path_rank.py`
  - root-to-tip style ranking from `marker_status.tsv`
  - weighted score terms:
    - derived support (`+1` default)
    - ancestral conflict (`-1` default)
    - nocall (`0` default)
    - derived `C>T` / `G>A` optional down-weight (`0.5` default)

Example commands:
```bash
annotate/y_haplo_from_markers.py \
  -i <output_dir>/iceman.vcf \
  --markers /tmp/markers_small.vcf.gz \
  -o /tmp/iceman_markers_vcf2 \
  --site-filter-mode deepvariant \
  --min-gq 20 --min-dp 3

annotate/y_path_rank.py \
  --input /tmp/iceman_markers_vcf2.marker_status.tsv \
  --out /tmp/iceman_markers_vcf2.path_rank.tsv \
  --clade-prefix G
```

Interpretation note:
- A candidate that appears only in gVCF top-path output can be caused by one or a few `nocall -> ancestral` state changes from non-variant blocks. Treat low-score gVCF-only additions as weak ranking shifts unless supported by additional derived markers.

## Terminal Y Placement Run (2026-07-25)

Goal: resolve the terminal node rather than confirm a published label. Uses read-level pileup
on the final BAM as primary evidence; caller output only as cross-check.

### Inputs
- BAM: `<analysis_dir>/dedup_out50/iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam`
  (chrY: 5,642,643 mapped reads)
- Reference: `mapping/index/hg38p14DH3630O.fa`
- Markers: `resources/snps_hg38.vcf.gz` (YBrowse, 2,920,250 rows)
- Cross-check callsets: `<output_dir>/iceman.vcf` (small model on), `<output_dir>/iceman-nosmall.vcf.gz`

### Commands
```bash
# chrY slice of the callset
bcftools view -t chrY <output_dir>/iceman.vcf -Oz -o iceman.chrY.vcf.gz
bcftools index -t iceman.chrY.vcf.gz      # 33,546 chrY records

# marker states from the callset
annotate/y_haplo_from_markers.py -i iceman.chrY.vcf.gz \
  --markers resources/snps_hg38.vcf.gz -o iceman_vcf \
  --site-filter-mode deepvariant --min-gq 20 --min-dp 3

# G-branch marker positions -> reads -> pileup (54,823 positions)
awk -F'\t' '{gsub(/"/,"",$8); gsub(/"/,"",$7);
  if ($8 ~ /^G/ || $7 ~ /^G[-2]/) print $1"\t"($2-1)"\t"$2}' \
  iceman_vcf.marker_status.tsv | sort -k2,2n > gall.bed
samtools view -b -M -L gall.bed "$BAM" -o gall.reads.bam && samtools index gall.reads.bam
samtools mpileup -f "$REF" -l gall.pos -r chrY -q 25 -Q 20 -d 1000 --no-BAQ gall.reads.bam \
  > gall.mpileup.txt
```

Ops notes:
- `samtools mpileup -a` with `-r chrY` walks all 57 Mb regardless of `-l`; it does not subset. Drop `-a`.
- `bcftools mpileup` on this BAM splits into 36 samples by `@RG` `SM`. Pass `--ignore-RG` or use
  `samtools mpileup` when a single merged pileup is wanted.
- Extracting reads for the marker BED first (`samtools view -M -L`) makes the pileup ~12 s instead of minutes.

### Results
Panel-wide over 52,632 covered G-branch markers: `ancestral 51,569`, `DERIVED 450`, `mixed 436`,
`other 146`, `no_coverage 2,222`.

Node-defining SNP sets (YFull tree, checked 2026-07-25):

| node | SNPs tested | result |
|---|---|---|
| backbone | `P287/PF3140`, `L91/S285/PF3246`, `PF3239`, `Z6043` | all derived, 0 ancestral reads |
| `G-L166` | `L166`, `L167`, `Z6516/FGC5675`, `FGC5696`, `FGC5721`, `Z6208`, `Z6219`, `Z6287`, `S19530/Z6213` | all derived, 0 ancestral reads (5 transversions) |
| `G-Z6494` | `Z6494/FGC5674` (DP 10), `FGC5687` (DP 3), `Z6215` (DP 4) | all ancestral, 0 derived reads |
| `G-Z6494` subclades | `Z6211` (DP 11), `Z6495/FGC5722` (DP 10), `FT84409` (DP 2) | all ancestral |

Caller cross-check at the five G-L166 SNPs present as variants — identical in both DeepVariant runs
(`PASS`, `GT=1/1`): `8525805`, `13776249`, `16280455`, `19013205`, `21843737`.

gVCF cross-check (`bcftools view -t chrY <output_dir>/iceman.gvcf`, 192,020 chrY records):
- G-L166 SNPs: same five `PASS` / `GT=1/1` rows.
- G-Z6494 branch SNPs (`FGC5687`, `Z6215`, `Z6211`, `Z6494`): all fall inside non-variant blocks with
  `GT=0/0`, `GQ=50`. DeepVariant calls them homozygous reference, not no-call — so "ancestral" here is a
  positive call, not absence of data.
- Note `MIN_DP` on those rows (2) is the block minimum over a 500-1000 bp span, not the depth at the
  marker. Site depth comes from the pileup (e.g. `Z6494` DP 10). Use the pileup for per-marker depth.

**Assignment: `G-L166*` (ISOGG-style `G2a2a1a2a1a*`).**
`G-Z6494` is the only child of `G-L166` in the current YFull tree, and the sample is ancestral for
all three of its defining SNPs. The published `G2a2a1a2a1a1b (G-Z6208*)` label follows ISOGG's
provisional (`~`) placement of `Z6208`, which the read data contradicts. See NOTES.md.

Saved artifacts:
- `results/iceman_y_L166_evidence.tsv`
- `results/iceman_y_deep_G_subtree.tsv`

## FTDNA Block Tree Walk (2026-07-25)

Purpose: test blocks copied from FTDNA's Block Tree at read level, since ISOGG does not rank most of
these markers and the label-driven path could not reach them.

One-time index build (~2 min, 796 MB GFF3 -> 3,121,999 rows):
```bash
annotate/build_marker_index.sh          # -> resources/marker_index.tsv.gz
```

Per-block test (~6 s for 24 markers):
```bash
BAM=<analysis_dir>/dedup_out50/iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam
REF=mapping/index/hg38p14DH3630O.fa
annotate/y_markers_pileup.py --bam "$BAM" --ref "$REF" \
  --marker-file block_pf3239.txt --label G-PF3239 --out results/block.tsv
```
Defaults: `-q 25 -Q 20 -d 1000 --no-BAQ`, 300 bp read-extraction flank. Markers absent from the
catalogue are reported as `not_in_catalogue`, not dropped.

Blocks tested (all reproduce the earlier hand-built pileups exactly):

| block | markers | derived | ancestral | mixed | verdict |
|---|---|---|---|---|---|
| `G-Z6488` | 3 | 3 | 0 | 0 | on path |
| `G-PF3238` | 5 | 5 | 0 | 0 | on path |
| `G-Z6128` | 2 | 2 | 0 | 0 | on path |
| `G-FT156872` | 1 | 0 | 1 | 0 | **excluded** (sibling) |
| `G-PF3239` | 24 | 22 | 0 | 2 | on path |
| `G-FGC2315` | 1 | 0 | 1 | 0 | **excluded** (sibling) |

Two `mixed` calls in `G-PF3239` are paralogous, not real heterozygosity: `Z6309` chrY:11265625 at
9 anc / 8 der on DP 17, `Z6312` chrY:11436604 at 2/4. Haploid chromosome — treat balanced sites as
collapsed-repeat artifacts.

Also re-ran the marker script to pick up the new `unplaceable_derived` output (49 s):
```bash
annotate/y_haplo_from_markers.py -i iceman.chrY.vcf.gz --markers resources/snps_hg38.vcf.gz \
  -o recheck --site-filter-mode deepvariant --min-gq 20 --min-dp 3
# total_markers: 2920250 / derived: 1569 / derived_unplaceable: 125
```

Saved artifacts:
- `results/iceman_y_ftdna_block_evidence.tsv`
- `results/iceman_y_unlabelled_derived_markers.tsv`

Still upstream of `L166` (`PF3239` = `G2a2a1a2a1`, `L166` = `G2a2a1a2a1a`). Open: what `G-L166` splits
into on FTDNA's tree.

## Swiss Neolithic L166 Test (2026-07-25)

Pre-registered in `PREREG_swiss_neolithic_L166.md` before any read of this dataset was examined.
Tests whether the Oberbipp/Rapperswil Horgen-context males published as terminal `PF3239` are in
fact derived at `L166` — i.e. whether the "7 of 10 are Ötzi's clade" claim circulating online is
a real result or an aggregator relabelling one node too deep.

Staging (9 runs, 607 MB, MD5-verified against ENA):
```bash
annotate/fetch_ena_runs.sh PRJNA608699 /mnt/AncientDNA/SwissLN-2020 \
  '^(MX182|MX187|MX209|MX210|MX211|MX212|MX213|SX8|SX10)$'
```

> **Correction, 2026-07-26.** This command was recorded with accession `PRJEB36959`, which is
> **not** this study — ENA resolves it to *"Partial sequence of 16S Leptospira borgpetersenii and
> Leptospira interrogans"*, and it returns zero read runs. The staged files' accessions
> (`SRR111793xx`, preserved in `manifest.tsv`) belong to `PRJNA608699`, which is what the run
> actually used. The command as published would have failed outright rather than fetched the wrong
> data, so no result is affected — but it did not reproduce, which is the point of writing it
> down. Second instance of documentation drift in this file; see the DeepVariant filename note
> above. The staging directory's `manifest.tsv` is the authoritative record, not this block.

**Gate: tool validation before any Swiss read was mapped** (5.6 s, 21 sites):
```bash
BAM=<analysis_dir>/dedup_out50/iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam
annotate/y_sites_pileup.py --bam "$BAM" --ref mapping/index/hg38p14DH3630O.fa \
  --sites markers/iceman_novel_candidates_all21.tsv --sample Iceman \
  > results/iceman_y_novel_candidates_regen.tsv
```
Reproduces `results/iceman_y_novel_branch_candidates.tsv` at **21/21 verdicts**, REJECT/MARGINAL
classification included. Five cells differed; all five were errors in the older hand-built table
and have been corrected in it. See `NOTES.md` for the specifics. The prereg commits to this gate
being blocking: a failure to reproduce would have been the finding, and nothing downstream would
have run.

Mapping (sequential, `-t 8`; bwa aln peak RSS ~4 GB, so one job at a time):
```bash
REF=$PWD/mapping/index/hg38p14DH3630O.fa \
  mapping/map_se_batch.sh /mnt/AncientDNA/SwissLN-2020 8
```
56 bp single-end 1240k capture, so `bwa aln`/`samse` at nf-core/eager 2.5.0 defaults
(`-n 0.01 -k 2 -l 1024`) rather than `bwa mem`. Same reference as the Iceman run. Dedup via
`samtools markdup -r` is mandatory here — the decision rules count independent reads.

Genotyping, once BAMs exist (marker sets committed under `markers/`):
```bash
REF=$PWD/mapping/index/hg38p14DH3630O.fa \
  annotate/y_genotype_batch.sh /mnt/AncientDNA/SwissLN-2020/bam results/swiss
```
Runs every marker set and the novel-site list across every BAM, adding the sample column neither
pileup tool emits on its own, and collates `swiss_coverage.tsv` with chrY DoC on the 12.57 Mb
callable denominator.

`backbone_control.txt` carries the positive control: Furtwängler 2020 reports these individuals
derived at `PF3239`, so independent recovery validates the whole chain before any `L166` call is
believed.

**Driver validated against the Iceman BAM before any Swiss BAM existed.** Run through the same
`y_genotype_batch.sh` with the BAM symlinked as `Iceman.rmdup.bam`, it reproduces the established
result exactly:

| set | markers | result |
|---|---|---|
| `L166_defining` | 9 | **9 DERIVED**, 0 ancestral reads at any (5 transversions, DP 2-10) |
| `Z6494_exclusion` | 3 | **3 ancestral**, 0 derived reads at any (DP 3-10) |
| `backbone_control` | 10 | **10 DERIVED**, 0 ancestral reads at any |

That is the whole G-L166-terminal finding recovered end to end by committed tooling, so a Swiss
sample's calls can be read against a chain that is known to work rather than one assumed to.

Verdict:
```bash
annotate/y_prereg_verdict.py --dir results/swiss
```
Applies the registered H1/H2/H0 rule per sample. Non-informative calls cannot produce H2 — only
≥2 ancestral reads at a site can — so an uncovered sample lands on H0 by construction rather than
by the reader's restraint. Checked against both known answers: Iceman → `H1_derived`, MX182 →
`H0_no_power`.

**Running `y_markers_pileup.py` on the first Swiss sample exposed that it did not implement the
registered call rules at all** — no depth threshold, no damage rule, so it reported `L166` as
*ancestral* off a single DP-1 read. Fixed in `1010ac0`; see `NOTES.md`. No Iceman conclusion
moves. This is the reason the validation gate is run against a known answer before new data.

### Extension to 15 samples (2026-07-26) — `Z6219` localisation and Family A pooling

Registered first, in `PREREG_swiss_neolithic_L166.md` **Amendment 2**, sections C–E, written and
committed before any of the added BAMs existed.

Staging (6 runs, 132 MB, MD5-verified):
```bash
annotate/fetch_ena_runs.sh PRJNA608699 /mnt/AncientDNA/SwissLN-2020 \
  '^(MX150|MX183|MX204|MX219|MX299)$'
annotate/fetch_ena_runs.sh PRJNA608699 /mnt/AncientDNA/SwissLN-2020 '^MX203$'
```
`fetch_ena_runs.sh` was changed for this run: it previously wrote `manifest.tsv` with `>`, so
staging a second subset into the same directory would have **erased the provenance of the nine
files already there**. It now merges into the cumulative manifest, deduplicated on
`run_accession`, and downloads only the newly requested rows. `manifest.tsv` now carries all 15.

Reference table for who these are — from `haplogroup.info/all-ancient-dna.txt`, retrieved
2026-07-26 and kept as `/mnt/AncientDNA/all-ancient-dna.2026-07-26.txt` (the earlier session read
this file without saving a copy, so the queries behind the `NOTES.md` claims were not re-runnable;
the dated filename is the fix):

| sample | ISOGG2019 | YFull | autosomal | Family A? |
|---|---|---|---|---|
| `MX219` | `G2a2a` | `G-PF3147` | 0.048 | yes |
| `MX150` | `G2a2a1a2` | `G-L91` | 0.200 | yes |
| `MX204` | `G2a2a1a2a` | `G-Z6484` | 0.433 | **no (`n/a`)** |
| `MX183` | `G2a2a1a2a1` | `G-PF3239` | 0.808 | yes |
| `MX203` | *(none)* | *(none)* | 0.142 | yes (9th member) |
| `MX299` | *(none)* | *(none)* | 0.068 | **no (`n/a`)** |

Only `MX204` and `MX299` carry outgroup information; the other four are the same Y chromosome as
samples already genotyped. `MX204` at 0.433x autosomal and published `Z6488` is the one sample in
the study that can plausibly separate the two `Z6219` explanations.

Mapping (resumable; skips the nine existing BAMs):
```bash
REF=$PWD/mapping/index/hg38p14DH3630O.fa \
  mapping/map_se_batch.sh /mnt/AncientDNA/SwissLN-2020 8
```

### Unhedged `G-L166` set (2026-07-26)

Pre-registered in `PREREG_unhedged_L166.md`, committed `e37b8ab` **before any read was staged**.
Five individuals from three unrelated studies, all carrying the *unhedged* ISOGG longhand
`G2a2a1a2a1a` with bare `G-L166` in the YFull and Y-Haplotree columns.

The compilation returns seven such rows; two of them are one man. `UNTA58_68Sk1`
(MittnikScience2019) and `E09538` (OlaldeNature2018) share radiocarbon lab number **MAMS-29075**
(3870±30 BP), mtDNA `J1c` and coordinates, and the Mittnik ID decodes to the Olalde colloquial
description ("Feature 68 Skeleton 1" = site UNTA58, feature 68, skeleton 1). `SX10` is the
seventh and was already genotyped in `results/swiss15/`. **Any tally over that list is one high
unless the duplicate is collapsed.**

Panel membership was settled before staging, from the panel definition rather than from the data
(`results/panel/1240k_marker_membership.tsv`):

```bash
annotate/panel_membership.py \
  --panel /mnt/MyGenome/Genos/FuQ/51.2.2M.snp \
  --markers markers/L166_defining.txt markers/Z6494_exclusion.txt markers/backbone_control.txt \
  --chain /mnt/GenomicData/OpenSNP/puller/hg38ToHg19.over.chain \
  --out results/panel/1240k_marker_membership.tsv
```

18 of 22 markers on-panel. `L166`, `L167`, `Z6219` and `PF3239` are deliberate targets; `Z6494`,
`Z6208`, `Z6488` and `FGC5687` are not. Matching is by lifted position, not by name — name
matching alone returns 5 of 22 and is **wrong**, because the panel carries synonyms (`P15` as
`PF3112`, `M3308` as `L1259`, `L91` as `S285`, `P287` as `rs4116820`). All 18 positional hits
agree on both alleles, which is what validates the liftover. Panel md5 `6389a9c8…`.

Staging (10 runs, 1.3 GB, MD5-verified, three studies into one directory; the manifest merge
added in the previous run is what makes that safe):
```bash
D=/mnt/AncientDNA/Unhedged-2026
annotate/fetch_ena_runs.sh PRJEB23635 $D '^I5118$'
annotate/fetch_ena_runs.sh PRJEB35980 $D '^(I15942|I14677|I14678)$'
annotate/fetch_ena_runs.sh PRJEB34400 $D '^UNTA58_68Sk1$'
```

Study accessions were **verified by ENA query before use, not assumed**. Two plausible guesses
made while locating these resolved to a *Medicago truncatula* rhizosphere GWAS (`PRJEB25849`) and
a feline cardiomyopathy study (`PRJEB27187`). This is the same failure mode as the `PRJEB36959`
correction above.

| sample | study | runs | total bases | autosomal | site, date |
|---|---|---|---|---|---|
| `I5118` | `PRJEB23635` | `ERR2207344`, `ERR2207549` | 375.5 Mbp | 1.631 | Mezőcsát-Hörcsögös, HU; 3300–3000 BCE |
| `I14678` | `PRJEB35980` | `ERR3800865`, `ERR3800866` | 848.8 Mbp | 3.867 | Serra Crabiles, Sardinia |
| `I14677` | `PRJEB35980` | `ERR3800863`, `ERR3800864` | 717.2 Mbp | 3.204 | Serra Crabiles, Sardinia |
| `I15942` | `PRJEB35980` | `ERR3800873`, `ERR3800874` | 315.0 Mbp | 0.957 | Anghelu Ruju, Sardinia |
| `UNTA58_68Sk1` | `PRJEB34400` | `ERR3518170`, `ERR3518171` | 166.6 Mbp | 1.040 | Augsburg, Bavaria |

`I14677` and `I14678` are from the same tomb at identical coordinates and count as **one**
independent observation unless their reads show otherwise — see the prereg's Independence section.

Mapping:
```bash
REF=$PWD/mapping/index/hg38p14DH3630O.fa \
  mapping/map_se_batch.sh /mnt/AncientDNA/Unhedged-2026 8
```

**`mapping/map_se_batch.sh` and `mapping/map_se_adna.sh` were fixed before this run** (`09754a8`).
Every sample here has two runs, and the batch driver keyed its work off individual FASTQs while
naming output `<sample>.rmdup.bam` and skipping when that file existed — so the second run of each
sample would have been silently dropped, no error and no log line. Glob order made the loss uneven:
for `UNTA58_68Sk1` the smaller run sorts first, so 2.54M of 3.30M reads (77%) would have gone
missing from the shallowest sample in the set. Runs are now grouped by sample and merged **before**
`markdup`, because PCR duplicates of one library are duplicates whichever run they were sequenced
in. Regression: re-running the batch over `/mnt/AncientDNA/SwissLN-2020` reports 15 samples, 1 run
each, `mapped=0 skipped=15 failed=0` — those results are untouched.

## Open Questions
- Quantify survivor vs removed read characteristics beyond prefix counts (MAPQ distribution, context near repetitive loci).
- Decide whether alt/decoy calls are retained only as technical appendix or excluded from variant interpretation entirely.
- Decide whether gVCF is required for all exploratory runs (skip when not needed).

### Marker sensitivity / mappability sweep (2026-07-26)

Tool: `annotate/y_mappability.py` (commit `e95244e`, contig-name fix `1a2b`-series). Answers "can a
read at this marker come back at all", which no site-level QC can see: `ylib.mapq_audit()` reports
MAPQ at the site and is therefore blind to reads that left for another contig.

    ./annotate/y_mappability.py \
      --markers markers/L166_defining.txt markers/Z6494_exclusion.txt markers/backbone_control.txt \
      --source mapping/index/hg38p14DH3630O.fa \
      --target working=mapping/index/hg38p14DH3630O.fa \
      --target noalt=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \
      --target chm13=mapping/index/chm13v2.0_maskedY_rCRSDH3630O.fa \
      --target hs37d5=/mnt/GenomicData/hs37d5/bwa/hs37d5.fa \
      --read-lengths 35,45,60,90 --threads 4 \
      --out results/mappability/y_marker_mappability.tsv

Output: `results/mappability/y_marker_mappability.tsv` (352 rows), `run.log`. Runtime ~25 min under
load average ~10; four index loads, one per reference.

Reference notes recorded at the time:
- `chm13v2.0_maskedY_rCRSDH3630O.fa` carries T2T chrY at **62,460,029 bp** (GRCh38: 57,227,415) and
  zero `chrY_*` GRCh38 patch contigs. T2T-CHM13v2.0 is CHM13 — a hydatidiform mole with no Y — plus
  a finished chrY from **GIAB HG002/NA24385**, a different individual from GRCh38's BAC-derived Y and
  from a different haplogroup. A "reference allele" on that chrY is therefore one living man's
  genotype, which is why the `ref_matches_anc` column matters more for CHM13 than for GRCh38.
- `GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna` is the reference the repo's own alt-contig
  list is built against, so it is an exact control for the custom build rather than an approximation.

Findings in `NOTES.md` under "Sensitivity Test: What Fraction Of Reads Can Even Come Back".

### Allele-aware mappability, damage forensics, reproducibility (2026-07-26, later)

All findings are in `NOTES.md`; this section records how to regenerate them.

**Reference composition audit** — no oral decoy present; chrY fix-patch alignment blocks:

    awk '{print $1}' mapping/index/hg38p14DH3630O.fa.fai | sed -E 's/[0-9]+$//' | sort | uniq -c
    awk -F'\t' '$3=="chrY"{...}' mapping/additional_hg38_p14_A3630_contigs.alt   # 4 blocks
    for b in .../I14677.rmdup.bam ...; do samtools idxstats "$b"; done            # reads/kb per contig

**Deposit provenance** — every ENA "fastq" is generated from a submitted BAM:

    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERR14752008&result=read_run\
&fields=run_accession,read_count,submitted_format,fastq_ftp&format=tsv"
    samtools view -H https://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14752008/CGG017683.bam
    for b in /mnt/AncientDNA/*/bam/*.rmdup.bam; do samtools idxstats "$b"; done   # unmapped fraction

**1000G phase 3 chrY membership** (`--sites` added to `panel_membership.py` for this):

    zcat /mnt/GenomicData/1KG/1KGchrY/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz \
      | awk -F'\t' '!/^#/{printf "%s\t24\t0.0\t%s\t%s\t%s\n",($3=="."?"snp_24_"$2:$3),$2,$4,$5}' > 1kg_chrY.snp
    annotate/panel_membership.py --panel 1kg_chrY.snp \
      --sites markers/iceman_novel_candidates_all21.tsv \
      --chain /mnt/GenomicData/OpenSNP/puller/hg38ToHg19.over.chain

**Z6219 read forensics** (terminus positions and library damage rate):

    annotate/y_read_evidence.py --bam <bam> --ref mapping/index/hg38p14DH3630O.fa \
      --site chrY:13782251 --anc C --der T
    annotate/y_read_evidence.py --bam <bam> --damage-profile --region chrY --ends 8

**F6 recheck** — all L166 reads with no MAPQ floor, 15 Swiss BAMs:

    for b in /mnt/AncientDNA/SwissLN-2020/bam/*.rmdup.bam; do
      annotate/y_read_evidence.py --bam $b --ref mapping/index/hg38p14DH3630O.fa \
        --site chrY:21843737 --anc C --der A --min-mq 0
    done
    # 14 reads, 8 samples, all ancestral C, all MAPQ 37, all NM=0, zero derived at any MAPQ

**Allele-aware mappability** (F5, and the YFull L166 set):

    for al in anc der; do
      annotate/y_mappability.py --markers markers/yfull_L166_defining.txt \
        --source mapping/index/hg38p14DH3630O.fa \
        --target noalt=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \
        --read-lengths 45,100,150 --threads 4 --allele $al --out tb_$al.tsv
    done

**MAPQ mechanism diagnostic** — single reads, ancestral vs derived, showing `XT:A:U X0=1 X1>=1`:
build a 45 bp read centred on each marker, substitute the derived base, then

    bwa aln -n 0.01 -k 2 -l 1024 <ref> multi.fq | bwa samse <ref> - multi.fq

**Reproducibility regression** — defaults must reproduce committed tables byte-identically:

    REF=mapping/index/hg38p14DH3630O.fa PREFIX=swiss \
      annotate/y_genotype_batch.sh /mnt/AncientDNA/SwissLN-2020/bam <outdir> markers
    diff results/swiss15/swiss_*.tsv <outdir>/swiss_*.tsv    # all 5 IDENTICAL, 2026-07-26

### Not done, and explicitly outstanding

- ~~**chrY-only pipeline emulation — NEVER RAN.**~~ **Run 2026-07-26 and the prediction confirmed;
  see "chrY-only emulation" below.** The reason it had not run was also recorded wrongly here: the
  note said `bwa index` "did not complete under load", but the log stops after 2 minutes at 53% of
  the BWT and no `bwa` process survived, so it was **killed when its parent shell exited**, not slow.
  Relaunched under `setsid` it finished in **199 seconds**. Load was never the problem.
- **Test A (F1 check) — not run.** `Z6219` for `I5118` sits in `results/unhedged/unhedged_L166_defining.tsv`
  and has deliberately not been read; see `PREREG_Z6219_node.md` §9.
- **Test B (co-segregation scan) — not run.** Marker set committed as
  `markers/yfull_L166_defining.txt`; 15 usable positions of 32.
- **`Z6219`'s YFull rank (26/32) is not reproduced** by anything measured here: it is 1.000 at both
  alleles at 35/45/60/100/150 bp with zero MQ0 and zero off-target in every reference tested.
- **Uniqueness-based filtering is a proposal, not applied.** `--report-uniqueness` reports it;
  no call uses it. Changing the filter requires its own registration and a re-run of every dataset.

## Tests A and B of `PREREG_Z6219_node.md` (2026-07-26)

**Test A — F1 check.** A read of already-committed tables plus the §8 read-terminus evidence.

    awk -F'\t' 'NR==1 || $3=="Z6219" || $3=="L166"' results/unhedged/unhedged_L166_defining.tsv
    awk -F'\t' 'NR==1 || $3=="Z6219" || $3=="L166"' results/swiss15/swiss_L166_defining.tsv

Evidence file `results/z6219_node/I5118_read_evidence.txt` regenerated by:

    REF=mapping/index/hg38p14DH3630O.fa
    BAM=/mnt/AncientDNA/Unhedged-2026/bam/I5118.rmdup.bam
    annotate/y_read_evidence.py --bam $BAM --ref $REF --site chrY:13782251 --anc C --der T
    annotate/y_read_evidence.py --bam $BAM --ref $REF --site chrY:21843737 --anc C --der A
    annotate/y_read_evidence.py --bam $BAM --damage-profile --region chrY --ends 8

**Test B — co-segregation scan.** `markers/yfull_L166_defining.txt` was already in `markers/`, so
the committed batch script picks it up with no new arguments:

    REF=$PWD/mapping/index/hg38p14DH3630O.fa PREFIX=swiss \
      annotate/y_genotype_batch.sh /mnt/AncientDNA/SwissLN-2020/bam results/testB markers
    REF=$PWD/mapping/index/hg38p14DH3630O.fa PREFIX=unhedged \
      annotate/y_genotype_batch.sh /mnt/AncientDNA/Unhedged-2026/bam results/testB_unhedged markers

    annotate/y_pool_family.py --genotypes results/testB/swiss_yfull_L166_defining.tsv \
      --members markers/family_a_members.txt --pool-name FamilyA \
      > results/testB/swiss_yfull_pooled.tsv

    annotate/y_cosegregation.py \
      --calls results/testB/swiss_yfull_L166_defining.tsv \
      --calls results/testB/swiss_yfull_pooled.tsv \
      --calls results/testB_unhedged/unhedged_yfull_L166_defining.tsv \
      --upstream MX210,MX213,SX10,FamilyA --outgroup I14677,I14678 \
      --ref-marker Z6219 --markers markers/yfull_L166_defining.txt \
      --also I5118,I15942 --out results/z6219_node/cosegregation_yfull32.tsv

Result: `splits_with_Z6219=2` (`Z6219`, `FGC5671`), `stays_with_block=3`, `uninformative=17`.

**`FGC5671` follow-up** — damage forensics, duplicate check and allele-aware mappability:

    for s in MX210 MX213 SX10 MX209; do
      annotate/y_read_evidence.py --bam /mnt/AncientDNA/SwissLN-2020/bam/$s.rmdup.bam \
        --ref mapping/index/hg38p14DH3630O.fa --site chrY:7784648 --anc G --der A
      annotate/y_read_evidence.py --bam /mnt/AncientDNA/SwissLN-2020/bam/$s.rmdup.bam \
        --damage-profile --region chrY --ends 10
    done
    samtools view /mnt/AncientDNA/SwissLN-2020/bam/MX210.rmdup.bam chrY:7784648-7784648

    printf 'FGC5671\nZ6219\nL166\nFGC5672\n' > /tmp/fgc_markers.txt
    for al in anc der; do
      annotate/y_mappability.py --markers /tmp/fgc_markers.txt \
        --source mapping/index/hg38p14DH3630O.fa \
        --target noalt=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \
        --read-lengths 35,45,60 --threads 4 --allele $al --out fgc_$al.tsv
    done

`FGC5671` recovers 1.000 at both alleles at all three lengths, 0 MQ0, 0 off-target.

**F2 backbone check:**

    annotate/y_pool_family.py --genotypes results/testB/swiss_backbone_control.tsv \
      --members markers/family_a_members.txt --pool-name FamilyA \
      > results/z6219_node/familyA_backbone_pooled.tsv

**Reproducibility regression, both cohorts.** The Test B batch regenerated every previously
committed table; all ten are byte-identical:

    for f in L166_defining Z6494_exclusion backbone_control coverage novel_sites; do
      diff -q results/swiss15/swiss_$f.tsv     results/testB/swiss_$f.tsv
      diff -q results/unhedged/unhedged_$f.tsv results/testB_unhedged/unhedged_$f.tsv
    done

### Still outstanding after Tests A and B

- **`FGC5671` is not placed relative to `L166`.** It is derived in Oberbipp and ancestral in
  Sardinians, which does not distinguish "on the `Z6219` branch" from "a private Oberbipp variant".
  `I5118`, the only `L166`-derived man available, has **no coverage** at `FGC5671`. The claim that
  the node rests on two SNPs currently borrows YFull's block assignment for the second one.
- **F2 remains largely untested.** `MX210`, `MX213` and `SX10` have no `PF3239` coverage.
- ~~**The chrY-only pipeline emulation still has not run.**~~ Run 2026-07-26; see below.
- **Uniqueness-based filtering still not applied.** Unchanged; still needs its own registration.

## chrY-only pipeline emulation (2026-07-26) — prediction confirmed

Long outstanding, and the reason was misdiagnosed twice. `bwa index` was not slow: it was being
reaped with its parent shell. Detached properly it takes 199 s.

    SC=<scratchpad>
    setsid nohup bwa index $SC/chrYonly.fa > $SC/bwaidx.log 2>&1 < /dev/null &

    annotate/y_mappability.py \
      --markers markers/Z6494_exclusion.txt markers/L166_defining.txt \
      --source mapping/index/hg38p14DH3630O.fa \
      --target wholegenome=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \
      --target chrYonly=$SC/chrYonly.fa \
      --read-lengths 45 --threads 4 \
      --out results/mappability/chrYonly_emulation.tsv

Result table committed as `results/mappability/chrYonly_emulation.tsv`.

## E1 — `FGC5671` in the Iceman (2026-07-26)

Registered in `PREREG_Z6219_node.md` §Extensions before running; findings in `NOTES.md`.
One BAM only, fixed in the registration:

```bash
BAM=/mnt/mirrored/iceman-reanalysis/dedup_out50/iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam
REF=mapping/index/hg38p14DH3630O.fa

annotate/y_markers_pileup.py --bam "$BAM" --ref "$REF" \
  --marker-file markers/yfull_L166_defining.txt --label Iceman \
  --out results/z6219_node/iceman_yfull_L166_defining.tsv
```

Read-level and library evidence, required by §8 at any decisive C>T/G>A site. Note that
`--damage-profile` is a **mode switch**, not an addition — it suppresses the per-read dump, so the
two are separate invocations:

```bash
annotate/y_read_evidence.py --bam "$BAM" --ref "$REF" --damage-profile --region chrY --max-reads 200000
annotate/y_read_evidence.py --bam "$BAM" --ref "$REF" --site chrY:7784648  --anc G --der A
annotate/y_read_evidence.py --bam "$BAM" --ref "$REF" --site chrY:13321379 --anc C --der T
# > results/z6219_node/iceman_E1_read_evidence.txt
```

Allele-aware mappability, run at `Z6499` because an unmappable derived allele produces the same
picture as a genuinely ancestral man:

```bash
printf 'Z6499\nFGC5671\n' > /tmp/e1_markers.txt
for al in anc der; do
  annotate/y_mappability.py --markers /tmp/e1_markers.txt \
    --source mapping/index/hg38p14DH3630O.fa \
    --target noalt=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \
    --read-lengths 35,45,60,100 --threads 4 --allele $al --out results/z6219_node/e1_mappability_$al.tsv
done
```

The counting rule registered in advance, applied. Controls are the 9 markers of
`markers/L166_defining.txt`, already known derived in the Iceman and excluded as non-independent;
`site_qc pass` is required for a position to count either way, so the inherited 30% MQ0 threshold
is what makes the three high-MQ0 `DERIVED` calls untested rather than confirmed:

```bash
grep -v '^#' markers/L166_defining.txt | grep -v '^$' > /tmp/e1_control.txt
awk -F'\t' '
  NR==FNR{ctrl[$1]=1; next}
  FNR==1{for(i=1;i<=NF;i++)c[$i]=i; print "marker\tpos\tcall\tsite_qc\tE1_class"; next}
  {
    m=$c["marker"]; call=$c["call"]; qc=$c["site_qc"]; valid=(qc=="pass")
    if (m in ctrl) { klass="control(" (valid&&call=="DERIVED" ? "ok" : "FAIL") ")"; ctrlok+=(valid&&call=="DERIVED"); nctrl++ }
    else if (valid && call=="DERIVED")   { klass="block-confirmed"; conf++ }
    else if (valid && call=="ancestral") { klass="block-refuted";   ref++ }
    else                                 { klass="untested";        unt++ }
    printf "%s\t%s\t%s\t%s\t%s\n", m, $c["pos"], call, qc, klass
  }
  END{ printf "#\n# control: %d/%d DERIVED at site_qc pass (run valid iff 9/9)\n", ctrlok, nctrl
       printf "# block-confirmed=%d  block-refuted=%d  untested=%d  (of %d non-control)\n", conf, ref, unt, conf+ref+unt }
' /tmp/e1_control.txt results/z6219_node/iceman_yfull_L166_defining.tsv \
  > results/z6219_node/iceman_E1_classification.tsv
```

Result: `FGC5671` **DERIVED** 7/0 at 0% MQ0 — F3 does not obtain. Controls 9/9. Non-control
block-confirmed 7, **block-refuted 1 (`Z6499`, ancestral 10/0)**, untested 5.
