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

## Open Questions
- Quantify survivor vs removed read characteristics beyond prefix counts (MAPQ distribution, context near repetitive loci).
- Decide whether alt/decoy calls are retained only as technical appendix or excluded from variant interpretation entirely.
- Decide whether gVCF is required for all exploratory runs (skip when not needed).
