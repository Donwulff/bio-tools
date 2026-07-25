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

## Open Questions
- Quantify survivor vs removed read characteristics beyond prefix counts (MAPQ distribution, context near repetitive loci).
- Decide whether alt/decoy calls are retained only as technical appendix or excluded from variant interpretation entirely.
- Decide whether gVCF is required for all exploratory runs (skip when not needed).
