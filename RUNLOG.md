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
- `iceman-y/results/iceman_smallmodel_compare.summary.tsv`
- `iceman-y/results/iceman_smallmodel_compare.y_markers.tsv`

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
- `iceman-y/results/iceman_y_L166_evidence.tsv`
- `iceman-y/results/iceman_y_deep_G_subtree.tsv`

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
- `iceman-y/results/iceman_y_ftdna_block_evidence.tsv`
- `iceman-y/results/iceman_y_unlabelled_derived_markers.tsv`

Still upstream of `L166` (`PF3239` = `G2a2a1a2a1`, `L166` = `G2a2a1a2a1a`). Open: what `G-L166` splits
into on FTDNA's tree.

## Swiss Neolithic L166 Test (2026-07-25)

Pre-registered in `iceman-y/prereg/swiss_neolithic_L166.md` before any read of this dataset was examined.
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
  > iceman-y/results/iceman_y_novel_candidates_regen.tsv
```
Reproduces `iceman-y/results/iceman_y_novel_branch_candidates.tsv` at **21/21 verdicts**, REJECT/MARGINAL
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
  annotate/y_genotype_batch.sh /mnt/AncientDNA/SwissLN-2020/bam iceman-y/results/swiss
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
annotate/y_prereg_verdict.py --dir iceman-y/results/swiss
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

Registered first, in `iceman-y/prereg/swiss_neolithic_L166.md` **Amendment 2**, sections C–E, written and
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

Pre-registered in `iceman-y/prereg/unhedged_L166.md`, committed `e37b8ab` **before any read was staged**.
Five individuals from three unrelated studies, all carrying the *unhedged* ISOGG longhand
`G2a2a1a2a1a` with bare `G-L166` in the YFull and Y-Haplotree columns.

The compilation returns seven such rows; two of them are one man. `UNTA58_68Sk1`
(MittnikScience2019) and `E09538` (OlaldeNature2018) share radiocarbon lab number **MAMS-29075**
(3870±30 BP), mtDNA `J1c` and coordinates, and the Mittnik ID decodes to the Olalde colloquial
description ("Feature 68 Skeleton 1" = site UNTA58, feature 68, skeleton 1). `SX10` is the
seventh and was already genotyped in `iceman-y/results/swiss15/`. **Any tally over that list is one high
unless the duplicate is collapsed.**

Panel membership was settled before staging, from the panel definition rather than from the data
(`iceman-y/results/panel/1240k_marker_membership.tsv`):

```bash
annotate/panel_membership.py \
  --panel /mnt/MyGenome/Genos/FuQ/51.2.2M.snp \
  --markers markers/L166_defining.txt markers/Z6494_exclusion.txt markers/backbone_control.txt \
  --chain /mnt/GenomicData/OpenSNP/puller/hg38ToHg19.over.chain \
  --out iceman-y/results/panel/1240k_marker_membership.tsv
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
      --out iceman-y/results/mappability/y_marker_mappability.tsv

Output: `iceman-y/results/mappability/y_marker_mappability.tsv` (352 rows), `run.log`. Runtime ~25 min under
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
    diff iceman-y/results/swiss15/swiss_*.tsv <outdir>/swiss_*.tsv    # all 5 IDENTICAL, 2026-07-26

### Not done, and explicitly outstanding

- ~~**chrY-only pipeline emulation — NEVER RAN.**~~ **Run 2026-07-26 and the prediction confirmed;
  see "chrY-only emulation" below.** The reason it had not run was also recorded wrongly here: the
  note said `bwa index` "did not complete under load", but the log stops after 2 minutes at 53% of
  the BWT and no `bwa` process survived, so it was **killed when its parent shell exited**, not slow.
  Relaunched under `setsid` it finished in **199 seconds**. Load was never the problem.
- ~~**Test A (F1 check) — not run.**~~ **Run 2026-07-26; F1 does not obtain.**
- ~~**Test B (co-segregation scan) — not run.**~~ **Run 2026-07-26; F3 does not obtain**, on one
  marker (`FGC5671`), whose status as a shared branch marker was settled by E1.
- ~~**The registered outcome of `iceman-y/prereg/Z6219_node.md` has never been recorded.**~~ **Recorded
  2026-07-26 as `H1 supported`**; see §Registered outcome in that document.
- **`Z6219`'s YFull rank (26/32) is not reproduced** by anything measured here: it is 1.000 at both
  alleles at 35/45/60/100/150 bp with zero MQ0 and zero off-target in every reference tested, and
  E3 confirmed this holds under a uniqueness criterion as well.
- ~~**Uniqueness-based filtering is a proposal, not applied.**~~ **Tested and rejected 2026-07-26**
  (`iceman-y/prereg/uniqueness_filter.md`): it fails U1 and U2. `MAPQ >= 25` stands and the lower-bound
  caveat on every `L166`/`Z6208` derived count is now permanent.
- **F2 is untested, not passed.** `MX210`, `MX213` and `SX10` have no `PF3239` coverage. Closing it
  requires either deeper data at `PF3239` in those three or an explicit statement that it cannot be
  closed with the libraries in hand.
- ~~**Test C (21 hedged Aesch/Muttenz individuals) — not run**, and out of scope until it has its
  own power statement.~~ **Power statement written 2026-07-26**, `iceman-y/prereg/testC_aesch_muttenz.md`.
  The cohort is **15 individuals, not 21** — the other 6 hedged `G-L166*` are the Oberbipp men
  already in `iceman-y/results/swiss15/`. **No read has been staged.** Registered expectation: ~4.0 callable
  at `Z6219`, ~2.1 callable at both `Z6219` and `L166`, split power 0.79 at best and 0.23 at the
  lower bound of the marker-rate interval. Still the only registered test bearing on the
  single-site limitation of the current verdict.
- **Aesch has Oberbipp's kinship problem.** The 13 Aesch candidates fall into 4 documented families
  plus one unassigned man; `Aes1`, `Aes12`, `Aes19`, `Aes20`, `Aes21` and `Aes23` are all Family D.
  Seven kin groups across two sites is what Test C actually buys, not fifteen individuals.
- **Test C's cohort is not the set P3 named.** P3 says "published-`PF3239` males"; the registered
  cohort is defined by the compilation's re-derived `G-L166*` label. Furtwängler 2020 Supplementary
  Table 1 has **not** been re-read for Aesch, so the overlap is unknown. Listed as a required
  pre-staging check in `iceman-y/prereg/testC_aesch_muttenz.md` §7.

## Tests A and B of `iceman-y/prereg/Z6219_node.md` (2026-07-26)

**Test A — F1 check.** A read of already-committed tables plus the §8 read-terminus evidence.

    awk -F'\t' 'NR==1 || $3=="Z6219" || $3=="L166"' iceman-y/results/unhedged/unhedged_L166_defining.tsv
    awk -F'\t' 'NR==1 || $3=="Z6219" || $3=="L166"' iceman-y/results/swiss15/swiss_L166_defining.tsv

Evidence file `iceman-y/results/z6219_node/I5118_read_evidence.txt` regenerated by:

    REF=mapping/index/hg38p14DH3630O.fa
    BAM=/mnt/AncientDNA/Unhedged-2026/bam/I5118.rmdup.bam
    annotate/y_read_evidence.py --bam $BAM --ref $REF --site chrY:13782251 --anc C --der T
    annotate/y_read_evidence.py --bam $BAM --ref $REF --site chrY:21843737 --anc C --der A
    annotate/y_read_evidence.py --bam $BAM --damage-profile --region chrY --ends 8

**Test B — co-segregation scan.** `markers/yfull_L166_defining.txt` was already in `markers/`, so
the committed batch script picks it up with no new arguments:

    REF=$PWD/mapping/index/hg38p14DH3630O.fa PREFIX=swiss \
      annotate/y_genotype_batch.sh /mnt/AncientDNA/SwissLN-2020/bam iceman-y/results/testB markers
    REF=$PWD/mapping/index/hg38p14DH3630O.fa PREFIX=unhedged \
      annotate/y_genotype_batch.sh /mnt/AncientDNA/Unhedged-2026/bam iceman-y/results/testB_unhedged markers

    annotate/y_pool_family.py --genotypes iceman-y/results/testB/swiss_yfull_L166_defining.tsv \
      --members markers/family_a_members.txt --pool-name FamilyA \
      > iceman-y/results/testB/swiss_yfull_pooled.tsv

    annotate/y_cosegregation.py \
      --calls iceman-y/results/testB/swiss_yfull_L166_defining.tsv \
      --calls iceman-y/results/testB/swiss_yfull_pooled.tsv \
      --calls iceman-y/results/testB_unhedged/unhedged_yfull_L166_defining.tsv \
      --upstream MX210,MX213,SX10,FamilyA --outgroup I14677,I14678 \
      --ref-marker Z6219 --markers markers/yfull_L166_defining.txt \
      --also I5118,I15942 --out iceman-y/results/z6219_node/cosegregation_yfull32.tsv

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

    annotate/y_pool_family.py --genotypes iceman-y/results/testB/swiss_backbone_control.tsv \
      --members markers/family_a_members.txt --pool-name FamilyA \
      > iceman-y/results/z6219_node/familyA_backbone_pooled.tsv

**Reproducibility regression, both cohorts.** The Test B batch regenerated every previously
committed table; all ten are byte-identical:

    for f in L166_defining Z6494_exclusion backbone_control coverage novel_sites; do
      diff -q iceman-y/results/swiss15/swiss_$f.tsv     iceman-y/results/testB/swiss_$f.tsv
      diff -q iceman-y/results/unhedged/unhedged_$f.tsv iceman-y/results/testB_unhedged/unhedged_$f.tsv
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
      --out iceman-y/results/mappability/chrYonly_emulation.tsv

Result table committed as `iceman-y/results/mappability/chrYonly_emulation.tsv`.

## E1 — `FGC5671` in the Iceman (2026-07-26)

Registered in `iceman-y/prereg/Z6219_node.md` §Extensions before running; findings in `NOTES.md`.
One BAM only, fixed in the registration:

```bash
BAM=/mnt/mirrored/iceman-reanalysis/dedup_out50/iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam
REF=mapping/index/hg38p14DH3630O.fa

annotate/y_markers_pileup.py --bam "$BAM" --ref "$REF" \
  --marker-file markers/yfull_L166_defining.txt --label Iceman \
  --out iceman-y/results/z6219_node/iceman_yfull_L166_defining.tsv
```

Read-level and library evidence, required by §8 at any decisive C>T/G>A site. Note that
`--damage-profile` is a **mode switch**, not an addition — it suppresses the per-read dump, so the
two are separate invocations:

```bash
annotate/y_read_evidence.py --bam "$BAM" --ref "$REF" --damage-profile --region chrY --max-reads 200000
annotate/y_read_evidence.py --bam "$BAM" --ref "$REF" --site chrY:7784648  --anc G --der A
annotate/y_read_evidence.py --bam "$BAM" --ref "$REF" --site chrY:13321379 --anc C --der T
# > iceman-y/results/z6219_node/iceman_E1_read_evidence.txt
```

Allele-aware mappability, run at `Z6499` because an unmappable derived allele produces the same
picture as a genuinely ancestral man:

```bash
printf 'Z6499\nFGC5671\n' > /tmp/e1_markers.txt
for al in anc der; do
  annotate/y_mappability.py --markers /tmp/e1_markers.txt \
    --source mapping/index/hg38p14DH3630O.fa \
    --target noalt=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \
    --read-lengths 35,45,60,100 --threads 4 --allele $al --out iceman-y/results/z6219_node/e1_mappability_$al.tsv
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
' /tmp/e1_control.txt iceman-y/results/z6219_node/iceman_yfull_L166_defining.tsv \
  > iceman-y/results/z6219_node/iceman_E1_classification.tsv
```

Result: `FGC5671` **DERIVED** 7/0 at 0% MQ0 — F3 does not obtain. Controls 9/9. Non-control
block-confirmed 7, **block-refuted 1 (`Z6499`, ancestral 10/0)**, untested 5.

## E2 — `CGG017683` (Crimea) at the YFull `L166` positions (2026-07-26)

Registered in `iceman-y/prereg/Z6219_node.md` §Extensions before running; result H0, findings in `NOTES.md`.
Stage 1 only — a remote indexed fetch against the deposited GRCh37 BAM. Stage 2 (full re-map) was
authorised only on >= 2 reads at a decisive position and is therefore **not** run.

The deposited BAM's index is fetched once and kept outside the repository (it is 6.6 MB of binary
and `.gitignore` excludes `*.bai`; it briefly sat in the repo root by accident and was moved):

```bash
SC=<scratchdir>/cgg017683
curl -o $SC/CGG017683.bam.bai \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14752008/CGG017683.bam.bai
```

Lift the hg38 marker set down to GRCh37. The two `FT*` markers land in no same-strand chain block
and are excluded rather than given a coordinate; the remaining 20 are checked against the hs37d5
reference base before they reach a pileup:

```bash
annotate/y_lift_markers.py \
  --markers markers/yfull_L166_defining.txt \
  --index resources/marker_index.tsv.gz \
  --chain /mnt/GenomicData/OpenSNP/puller/hg38ToHg19.over.chain \
  --target-ref /mnt/GenomicData/hs37d5/hs37d5.fa --target-chrom Y \
  --sites-out $SC/sites_b37.tsv \
  --report-out iceman-y/results/z6219_node/e2_cgg017683_lift_report.tsv
# 22 markers: 20 lifted and ref-checked, 2 excluded (FT91632, FT191098: unmapped)
```

Fetch the chrY slice by remote index — 17 s, no bulk download — and genotype:

```bash
samtools view -b -X https://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14752008/CGG017683.bam \
  $SC/CGG017683.bam.bai Y -o $SC/chrY.b37.bam
samtools index $SC/chrY.b37.bam
samtools idxstats $SC/chrY.b37.bam   # Y 59373566 54685 0

annotate/y_sites_pileup.py --bam $SC/chrY.b37.bam --ref /mnt/GenomicData/hs37d5/hs37d5.fa \
  --sites $SC/sites_b37.tsv --sample CGG017683 --min-mq 25 --min-bq 20 \
  > iceman-y/results/z6219_node/e2_cgg017683_pileup.tsv
```

The zeros at the decisive positions are absence of molecules, not filtering. Checked with every
filter off, because this deposit ran `samtools calmd -Erb` without `-A` and a BAQ-capped base
quality would mimic no coverage:

```bash
for p in 15894131 23989884 23989903 15433259; do   # Z6219, L166, L167, Z6499 in b37
  samtools view -c $SC/chrY.b37.bam Y:$p-$p; done   # 0 0 0 0
samtools depth -a -q 0 -Q 0 -r Y:23989884-23989903 $SC/chrY.b37.bam   # 0
```

Result: 3 of 20 covered (15.0% vs 17.5% expected at 0.192x), all single-read. `FGC5721` derived
(registered `DERIVED_1read_transversion`, but uninformative — the Iceman is derived there too),
`Z6134` and `S19530` single-read nocalls. `Z6219`/`L166`/`L167`/`Z6499` zero reads. **H0.**

## E3 — uniqueness filter vs `MAPQ >= 25` (2026-07-26)

Registered in `iceman-y/prereg/uniqueness_filter.md` before running; **rejected**, findings in `NOTES.md`.

`y_mappability.py` gained two **reported-only** columns, `n_unique` and `frac_recovered_unique`,
counting what `XT:A:U && X0 == 1` would have kept. They sit behind the same correct-contig and
exact-CIGAR gates as `n_mq_ge_min`, so the two criteria share a denominator. Existing columns are
computed exactly as before; no call changes by running this. `iceman-y/results/mappability/y_marker_mappability.tsv`
predates the new columns and was **not** regenerated — its values are unaffected and no claim needs it.

```bash
for al in anc der; do
  annotate/y_mappability.py \
    --markers markers/L166_defining.txt markers/Z6494_exclusion.txt markers/backbone_control.txt \
    --source mapping/index/hg38p14DH3630O.fa \
    --target working=mapping/index/hg38p14DH3630O.fa \
    --target noalt=mapping/index/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna \
    --read-lengths 35,45,60,100 --threads 4 --allele $al \
    --out iceman-y/results/mappability/e3_uniqueness_$al.tsv
done
```

U1 and U2 applied exactly as registered — thresholds 0.05 and "must not exceed", both fixed before
the run:

```bash
awk -F'\t' '
 NR==FNR{ if(FNR==1){for(i=1;i<=NF;i++)c[$i]=i; next}
          k=$c["marker"]"|"$c["read_len"]"|"$c["target"]
          A_mq[k]=$c["frac_recovered"]; A_un[k]=$c["frac_recovered_unique"]
          A_off[k]=$c["n_offtarget"]; next }
 FNR==1{next}
 { k=$c["marker"]"|"$c["read_len"]"|"$c["target"]
   du=A_un[k]-$c["frac_recovered_unique"]; if(du<0)du=-du
   if(du>0.05) u1++                                                  # U1
   if($c["frac_recovered_unique"]<$c["frac_recovered"]) u1b++        # U1, derived worse
   if(A_off[k]+0>0 && A_un[k]>A_mq[k]) u2++                          # U2, ancestral tile
   if($c["n_offtarget"]+0>0 && $c["frac_recovered_unique"]>$c["frac_recovered"]) u2++
 } END{print "U1 fails:",u1+0,"  U1 derived-worse:",u1b+0,"  U2 fails:",u2+0}
' iceman-y/results/mappability/e3_uniqueness_anc.tsv iceman-y/results/mappability/e3_uniqueness_der.tsv
```

Result over 176 marker/length/target cells: **U1 fails 18** (28 were asymmetric under MAPQ, so the
bias is reduced, not repaired; `L166` at 45 bp goes 0.156 → 0.689 against a 0.05 threshold),
**U1 derived-worse 0**, **U2 fails 29** — including `FGC5687` at 60 bp (0.000 → 0.250) and 100 bp
(0.160 → 0.890) with 23 and 5–6 off-target tile reads. Not adopted. U3 not run: the decision was
already determined by U1 and U2.

## Test C power statement (2026-07-26) — no read staged

Pre-registered in `iceman-y/prereg/testC_aesch_muttenz.md`, committed together with the tables below.
**Nothing here touches Aesch or Muttenz sequence data**; the whole estimate is built from read
counts, an external coverage proxy, and the marker rates already measured on the Oberbipp cohort.

Cohort enumeration from the All Ancient DNA compilation. The file is latin1 with bytes that make
`grep` treat it as binary, so `-a` is required or every count silently comes back empty:

    F=/mnt/AncientDNA/all-ancient-dna.2026-07-26.txt
    grep -ac 'L166' "$F"          # 37

    awk -F'\t' -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}
      $c["Location"] ~ /Aesch|Muttenz/ && $c["Y-Haplotree-Variant"]=="G-L166*" {
        id=$c["Object-ID"]; ena=id; sub(/^Aes/,"Aesch",ena); sub(/^RA/,"SNPRA",ena)
        print ena,$c["NRY"],"compilation_id="id";site="$c["Location"]";published_Y="$c["Y-DNA"] }' \
      "$F" | sort -t$'\t' -k1,1V > iceman-y/results/testC_power/candidates.tsv

Calibration proxy (same columns, `Location ~ /Oberbipp|Rapperswil/`) →
`iceman-y/results/testC_power/proxy_calibration.tsv`. `iceman-y/results/testC_power/candidates_lineage.tsv` is hand-
built from the compilation's `Kinship-Notes` and applies the §8 pooling rule: only the documented
1st-degree male–male pair `Aes12`–`Aes19` is merged.

ENA availability and protocol match, confirmed for all 15 candidates and all 97 runs:

    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA608699\
&result=read_run&fields=run_accession,sample_alias,read_count,library_strategy,instrument_model&format=tsv"

All 15 are present, all `Targeted-Capture` / `Hybrid Selection` on HiSeq 3000 — the same protocol
and instrument as the calibration cohort, which is what licenses transferring the marker rates.

Power:

    annotate/y_power_estimate.py \
      --calib-genotypes iceman-y/results/testB/swiss_yfull_L166_defining.tsv \
      --calib-coverage  iceman-y/results/swiss15/swiss_coverage.tsv \
      --candidates      iceman-y/results/testC_power/candidates_lineage.tsv \
      --proxy-calibration iceman-y/results/testC_power/proxy_calibration.tsv \
      --markers Z6219,L166 --out-prefix iceman-y/results/testC_power/testC_lineage

and once more over the unpooled 15 with all four markers that attract any reads, which is what
produces `iceman-y/results/testC_power/testC_{calibration,candidates,power}.tsv`:

    annotate/y_power_estimate.py \
      --calib-genotypes iceman-y/results/testB/swiss_yfull_L166_defining.tsv \
      --calib-coverage  iceman-y/results/swiss15/swiss_coverage.tsv \
      --candidates      iceman-y/results/testC_power/candidates.tsv \
      --proxy-calibration iceman-y/results/testC_power/proxy_calibration.tsv \
      --markers Z6219,L166,FGC5671,L167 --out-prefix iceman-y/results/testC_power/testC

Self-check — the same command with `--candidates` pointing at the calibration cohort predicts
**2.08** callable at `Z6219` (observed **2**) and **3.15** at `L166` (observed **4**).

Marker-rate interval bracket, which dominates every other uncertainty here:

    for s in 0.479 1.0 1.838; do ... --rate-scale $s ... ; done

giving expected callable at `Z6219` of **1.3 / 4.0 / 7.9** and split power at frequency 0.5 of
**0.23 / 0.79 / 0.98**.

**Gotcha recorded.** `y_power_estimate.py` first used a through-origin fit for the coverage proxy.
`MX182`'s compilation row carries `NRY = 2` against 9,148 mapped chrY reads — a defect noted in
`NOTES.md` on 2026-07-25 — and that one row alone raised the fitted ratio from 0.6514 to 0.7455, a
14% optimistic bias in every downstream power number. The tool now takes the median of per-sample
ratios and reports the through-origin value alongside it for comparison. No sample is dropped.

### Test C §7 pre-staging checks (2026-07-26) — still no read staged

    cd <scratch>
    curl -sL -o furtw_suppl.pdf "https://static-content.springer.com/esm/\
art%3A10.1038%2Fs41467-020-15560-x/MediaObjects/41467_2020_15560_MOESM1_ESM.pdf"
    pdftotext -layout furtw_suppl.pdf furtw_suppl.txt
    grep -n "Supplementary Table" furtw_suppl.txt     # Table 1 at l.220, Table 5 at l.454

25 pp., 631 lines of extracted text. Supplementary Table 1 (Y haplogroup for all male individuals)
and Table 5 (relatedness, first and second degree) are the two read. Outcome in
`iceman-y/prereg/testC_aesch_muttenz.md` §"§7 checks — completed" and `NOTES.md`: all three checks pass, the
cohort is a superset of the published-`PF3239` set, and the registered `Aes12`–`Aes19` pooling is
the only male–male first-degree pair the publication documents among the candidates.

The PDF is not committed — it is a publisher artefact reachable from the URL above, and the two
tables read from it are quoted where they are used.

### Test C staging and mapping (2026-07-26) — first reads of this cohort

Staged into the existing Swiss Late Neolithic directory rather than a new one: same study, same
project accession, and `fetch_ena_runs.sh` merges into `manifest.tsv` rather than overwriting it, so
the provenance of the 15 Oberbipp/Rapperswil runs already there is preserved. `manifest.tsv` now
carries 30 rows.

```bash
annotate/fetch_ena_runs.sh PRJNA608699 /mnt/AncientDNA/SwissLN-2020 \
  '^(Aesch1|Aesch4|Aesch6|Aesch7|Aesch12|Aesch13|Aesch14|Aesch17|Aesch19|Aesch20|Aesch21|Aesch22|Aesch23|SNPRA61|SNPRA62)$'
```

15 runs, 28,601,005 reads, 734.5 MB, all MD5-verified against the ENA-published sums
(`failures=0`). One run per sample, so the multi-run grouping defect fixed in `map_se_batch.sh`
(a sample sequenced across several runs having all but the first silently skipped) cannot apply to
this cohort.

```bash
REF=$PWD/mapping/index/hg38p14DH3630O.fa \
  mapping/map_se_batch.sh /mnt/AncientDNA/SwissLN-2020 8
```

The batch sees 30 samples and maps 15: the resume check skips every Oberbipp BAM already present,
so those are not rebuilt and their calls are untouched. Same reference and same `bwa aln`/`samse`
parameters as the Iceman and Swiss 15 BAMs, which is what makes the two cohorts comparable without
a liftover.

**Both commands were run detached (`setsid`) from a session runner rather than in the foreground.**
The staging download and the sequential mapping pass together exceed any single tool invocation's
timeout, and a fetch killed mid-file leaves a truncated FASTQ. That is recoverable — the MD5 resume
check refetches it — but it is recoverable only because the check exists, and it did happen once
here.

**Instrument check during mapping, no genotype read.** chrY yield of the first two finished BAMs,
against the registered prediction (`NRY` proxy times the fitted ratio 0.6514):

```bash
cd /mnt/AncientDNA/SwissLN-2020/bam
for s in MX210 SX10 Aesch12 Aesch13; do
  printf "%-8s " "$s"
  echo "$(samtools view -c -q 25 "$s.rmdup.bam" chrY) chrY@MQ25 / $(samtools view -c "$s.rmdup.bam" chrY) all"
done
```

| sample | proxy | predicted | observed @ MQ25 | all |
|---|---|---|---|---|
| `MX210` (fitted) | 38,316 | — | 30,629 | 31,737 |
| `SX10` (fitted) | 16,319 | — | 10,663 | 10,929 |
| `Aesch12` | 20,241 | 13,185 | 12,082 | 12,556 |
| `Aesch13` | 17 | 11 | **15,953** | 16,493 |

`Aesch12` is 8% under prediction, which validates the calibration out of sample. `Aesch13` is the
`MX182` defect repeating: a near-zero `NRY` against a normal library. See `NOTES.md`. §4 of the
pre-registration is **not** rewritten — the registered power stands as predicted, and is now known to
be conservative.

### Corroborating the `NRY` defect from the publication (2026-07-26) — still no genotype read

The machine-readable supplements were already in the scratchpad from an earlier session. `MOESM2` is
Supplementary Table 1 as a spreadsheet, with `SNPs 1240k`, `Mapped Reads after RMDup`, `sex`,
`1st degree relatives` and a `Y GH` column.

Read with a standalone parser rather than a dependency install (the box was out of swap at the time):

    cd <scratch>/swiss_supp
    python3 - <<'PY'
    import zipfile, re
    import xml.etree.ElementTree as ET
    NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
    def colnum(ref):                      # 'B3' -> 1. REQUIRED: xlsx omits empty
        n = 0                             # cells, so document order shifts rows left.
        for ch in re.match(r'([A-Z]+)', ref).group(1):
            n = n * 26 + (ord(ch) - 64)
        return n - 1
    z = zipfile.ZipFile('41467_2020_15560_MOESM2_ESM.xlsx')
    ss = [''.join(t.text or '' for t in si.iter(f'{NS}t'))
          for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall(f'{NS}si')]
    for r in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).iter(f'{NS}row'):
        row = {}
        for c in r.findall(f'{NS}c'):
            v = c.find(f'{NS}v')
            if v is not None and v.text is not None:
                row[colnum(c.get('r'))] = ss[int(v.text)] if c.get('t') == 's' else v.text
        print([row.get(i, '') for i in range(max(row) + 1 if row else 0)])
    PY

Columns of interest: 5 `ID`, 6 `1st degree relatives`, 13 `sex`, 15 `SNPs 1240k`, 23 `Y GH`.

Three results, all in `NOTES.md`: `Aes13` has **565,667** SNPs on the 1240k panel against a
compilation `NRY` of 17, which settles that defect; the spreadsheet's `Y GH` column contradicts the
PDF's Y-haplogroup table for eight of the fifteen candidates; and the `1st degree relatives` column
confirms `Aes12`–`Aes19` as the only within-cohort male–male first-degree pair, independently of
Supplementary Table 5.

`iceman-y/prereg/testC_aesch_muttenz.md` §7 check 3 was run against the PDF's Y table and **stands unamended**.
The `.xlsx` files are not committed, for the same reason the PDF is not: publisher artefacts,
reachable from the article's supplementary section, quoted where used.

### Test C §7 check 5 — backbone control, run alone before any `L166` table existed (2026-07-26)

Mapping finished `mapped=15 skipped=15 failed=0`. The 15 Test C BAMs were then exposed through a
symlink view so the Oberbipp BAMs in the same directory are not genotyped into Test C's tables:

```bash
cd /mnt/AncientDNA/SwissLN-2020 && mkdir -p bam-testC && cd bam-testC
for s in Aesch1 Aesch4 Aesch6 Aesch7 Aesch12 Aesch13 Aesch14 Aesch17 \
         Aesch19 Aesch20 Aesch21 Aesch22 Aesch23 SNPRA61 SNPRA62; do
  ln -sf "../bam/${s}.rmdup.bam" .; ln -sf "../bam/${s}.rmdup.bam.bai" .
done
```

The gate was run against a marker directory containing **only** `backbone_control.txt`, so that no
`L166`-defining call existed anywhere on disk at the moment the control was read:

```bash
G=<scratch>/markers_gate; mkdir -p "$G"; cp markers/backbone_control.txt "$G/"
REF=$PWD/mapping/index/hg38p14DH3630O.fa PREFIX=testC SITES= \
  annotate/y_genotype_batch.sh /mnt/AncientDNA/SwissLN-2020/bam-testC iceman-y/results/testC "$G"
```

150 cells (15 samples x 10 markers):

| call | n |
|---|---|
| `no_coverage` | 95 |
| `DERIVED_1read_transversion` | 23 |
| `DERIVED` | 17 |
| `nocall_damage_prone_1read` | 13 |
| `mixed` | 1 |
| `other_allele` | 1 |
| **`ancestral`** | **0** |

| marker | derived | ancestral | covered |
|---|---|---|---|
| `L91` | 10 | 0 | 10 |
| `M3308` | 9 | 0 | 9 |
| `P287` | 8 | 0 | 8 |
| `Z6043` | 7 | 0 | 8 |
| `P15` | 5 | 0 | 8 |
| `PF3147` | 1 | 0 | 3 |
| `PF3239` | **0** | 0 | 7 |
| `PF3148`, `PF3177` | 0 | 0 | 1 each |
| `Z6488` | 0 | 0 | 0 |

**Gate passes: 40 derived calls, zero ancestral calls anywhere.** The published path
`P15 > M3308 > ... > L91` is recovered independently in this cohort, and nothing contradicts it.

**But not via `PF3239`, which is what the marker file's own comment names as the positive control.**
`PF3239` is a C>T transition at chrY:15205748. Seven libraries cover it, every one of them with
**exactly one derived read and no ancestral read**; the other eight have no coverage. Under the
registered rules every one of those is `nocall_damage_prone_1read` and **is not evidence**. So the
chain is validated by `L91`, `M3308`, `P287`, `Z6043` and `P15`, and the specific marker the
publication assigns these men is uninformative here. That is a depth-rule consequence, not a
contradiction of the published call — 7 of 7 covered cells carrying the derived base with zero
ancestral reads is consistent with it, and is deliberately not counted as confirming it.

Two non-ancestral oddities, neither contradicting the chain: `Aesch21` at `PF3147` is `mixed`
(1 ancestral, 2 derived, C>T, `pct_mq0 = 0%`), and `Aesch4` at `Z6043` has a single read carrying
neither allele (`other_allele`, transversion, `pct_mq0 = 0%`).

### Test C genotyping (2026-07-27)

```bash
REF=$PWD/mapping/index/hg38p14DH3630O.fa PREFIX=testC \
  annotate/y_genotype_batch.sh /mnt/AncientDNA/SwissLN-2020/bam-testC iceman-y/results/testC
```

Four marker sets across 15 samples: `testC_L166_defining.tsv` (135 rows),
`testC_yfull_L166_defining.tsv` (330), `testC_backbone_control.tsv` (150),
`testC_Z6494_exclusion.tsv` (45), `testC_novel_sites.tsv` (120), plus coverage and params.

Damage evidence, required by §8 for the decisive C>T calls:

```bash
B=/mnt/AncientDNA/SwissLN-2020/bam-testC
R=/home/jsantala/src/bio-tools/mapping/index/hg38p14DH3630O.fa
for s in Aesch12 Aesch13 Aesch23 SNPRA61 SNPRA62; do
  annotate/y_read_evidence.py --bam "$B/${s}.rmdup.bam" --ref "$R" \
    --site chrY:13782251 --anc C --der T
  annotate/y_read_evidence.py --bam "$B/${s}.rmdup.bam" --ref "$R" \
    --damage-profile --region chrY --max-reads 200000
done
```

Outcome in `iceman-y/prereg/testC_aesch_muttenz.md` §"Registered outcome": PC1 met, FC1 and FC2 not triggered,
FC3 not evaluable, PC3 holds.

**Two gotchas, both mine, both recorded because each produced plausible output rather than an error.**

*The shell's working directory persists between tool calls.* A `cd iceman-y/results/testC` from an earlier
command was still in effect when a later block set `R=$PWD/mapping/index/...`, which expanded to a
path that does not exist. `samtools faidx` failed, `ref_base()` returned its `"?"` fallback, and the
run printed `ref=?` and silently dropped the damage-pattern tag from every verdict. The per-read
bases, positions and qualities come from the BAM and were unaffected, but the run was repeated with
an absolute reference path before anything was recorded from it. **Use absolute paths in these
commands; `$PWD` is not stable across calls.**

*Read-ID adjacency is not a duplicate signature here.* Every pair of derived reads at `Z6219` carried
consecutive SRA spot indices (`...1527837`/`...1527838` and so on) in all five samples, which looks
like undetected duplication. The control is to look at a site with more reads: `L166`'s ancestral
reads in the same libraries are consecutive too (`2687332`, `2687333`, `2687336`-`2687339`). ENA
regenerated these FASTQs from coordinate-sorted BAMs, so reads at one locus have adjacent spot
numbers by construction. The reads differ in length and strand and survive `samtools markdup -r`.

**`derived(damage-pattern)` in `y_read_evidence.py` output is mechanical**: it tags any C>T or G>A by
mutation class, not by how the read looks. It is not a damage assessment and carries no evidence on
its own; `dist_5p`/`dist_3p` against the library profile is what does.

## The tree gets a file — `markers/tree_local.tsv` and the placement tools (2026-07-27)

Topology had never been machine-readable. Marker infrastructure handled two things — named markers
resolved through the YBrowse index, and bare coordinates — and neither expresses which node sits
above which. That `Z6219` is above `L166` and `Z6499` below it existed only as prose in `NOTES.md`,
so no script could act on it and no script could be checked against it.

`markers/tree_local.tsv` now states it: one row per node, with a `status` column of
`published` / `provisional` / `putative` / `refuted`, a `defining` marker list (derived means this
node or below; a terminal call may rest only on these) and an `equivalent` list (block members no
sample here orders against the node — counted, reported, never decisive).

`annotate/ytree.py` loads and validates it, and carries `TreeScorer`, shared by both consumers for
the reason `ylib.py` exists. `annotate/y_tree_place.py` runs it against read-level pileup tables;
`y_path_rank.py --tree` runs it against VCF-derived marker_status tables alongside the existing
label ranking. The two are independent by construction — the ranking reads HG/ISOGG labels, the
scorer reads only marker names — so they can disagree, which is the point.

```bash
python3 annotate/test_ytree.py            # 40 checks, no pytest, no fixtures on disk

mkdir -p iceman-y/results/placement
python3 annotate/y_tree_place.py --sample Iceman \
  --pileup iceman-y/results/iceman_y_L166_evidence.tsv iceman-y/results/z6219_node/iceman_yfull_L166_defining.tsv \
  --out iceman-y/results/placement/iceman_placement.tsv
python3 annotate/y_tree_place.py --pileup iceman-y/results/swiss15/swiss_*.tsv \
  --out iceman-y/results/placement/oberbipp15_placement.tsv
python3 annotate/y_tree_place.py --pileup iceman-y/results/testC/testC_*.tsv \
  --out iceman-y/results/placement/aesch_muttenz_placement.tsv
python3 annotate/ytree.py --newick iceman-y/results/placement/tree_local.nwk \
  --markers-out iceman-y/results/placement/tree_local_markers.tsv
```

**The ladder now comes out of a tool instead of a paragraph.** The Iceman is `G-L166`, on two
derived transversions, with `G-Z6494`, `G-Z6211` and `G-Z6499` all ancestral below him; `G-Z6208` is
`DERIVED` and `refuted`, which is the published label's error in one line. Oberbipp gives `MX210`
and `SX10` at `G-Z6219` and nobody at `G-L166`; Aesch/Muttenz gives `Aesch12`, `Aesch13`, `Aesch23`,
`SNPRA61` and `SNPRA62` at `G-Z6219`. Those reproduce the registered outcomes of the Swiss run and
Test C, from the committed tables, without a human in the loop.

**Two defects the first live run exposed, both fixed before the tables above were written.**

*Overlapping block files double-counted.* `markers/L166_defining.txt` and
`markers/yfull_L166_defining.txt` share ten positions, so globbing a results directory fed the same
read in twice and every derived and ancestral total came out inflated. Rows are now deduplicated on
`(sample, marker)`, and a call that differs between two input files is reported rather than silently
resolved.

*`Aesch7` places at `G-L166`, and the tool has to say why that is thin.* He has one derived read at
`L167` (`DERIVED_1read_transversion`, a derived call under the standing rules) and one ancestral read
at `L166` (`low_power_1read_ancestral`, a nocall under the same rules). Both rules are correct and
together they produce a `G-L166` placement from one read. Test C reported the cohort aggregate — "24
ancestral against 1 derived" at `L167` — so this never surfaced per-sample. The scorer does not
change either rule; it counts low-power calls separately and prints `weak_ancestral=1` against the
placement. `Aesch7` is the only sample in either cohort carrying that caveat.

**`G-Z6208` is the reason the `refuted` status exists.** The published Iceman label
`G2a2a1a2a1a1b (G-Z6208*)` treats `Z6208` as a node below `L166`. The marker is real and the Iceman
is derived at it; the node is not. A `refuted` row keeps both facts — the call is scored and printed,
and the node can never appear on a path, be a terminal placement, or reach the Newick export.

## The investigation gets a directory, and a verdict page (2026-07-27)

The Iceman work had no home. Five pre-registrations, a protocol doc and eleven results directories
sat in the repository root beside `README.md`, `AGENTS.md` and `TODO.md`, which are tooling. That is
why a findings page had no good name: there was no context for it to sit in.

`iceman-y/` now holds the case-specific material — `prereg/` (the five registrations, `PREREG_`
prefix dropped since the directory says it), `PROTOCOL_extending_analyses.md`, and `results/` with
every table the verdict cites. Root `results/` keeps only the three QC PDFs, which are mapping-side.
`iceman-y/README.md` is the verdict, and needs no generic filename because the directory names it.

`NOTES.md`, `METHODS.md` and `RUNLOG.md` deliberately stayed at root. `AGENTS.md` describes them as
remapping methods and pairing/EAGER gotchas; they predate this investigation and are shared with the
mapping work. `NOTES.md` carries no `##` heading until line 1244. Splitting them would mean choosing
a boundary that is not there, across 110 cross-references, and the loss would be the mapping notes'
continuity rather than the Iceman notes' findability.

The move was `git mv` (history preserved by rename detection) plus a scripted rewrite of every path
reference across `*.md`, `*.py` and `*.sh` outside `DeDup/`. Root-relative paths were used
throughout rather than paths relative to each document, because every one of these references is a
backticked path or a shell command run from the repository root — there were **zero** markdown links
among them, checked before rewriting. 18 files, 199 lines.

Verified afterwards: no stale reference to any moved path remains; `annotate/test_ytree.py` passes
40/40; `y_tree_place.py` reproduces the Aesch/Muttenz placements from the new locations.

**One new table, and it closes the ladder.** The unhedged cohort had never been run through the
placement tool:

```bash
python3 annotate/y_tree_place.py \
  --pileup iceman-y/results/unhedged/unhedged_*.tsv \
           iceman-y/results/testB_unhedged/unhedged_*.tsv \
  --out iceman-y/results/placement/unhedged_placement.tsv
```

`I5118` places at `G-L166` — the only individual other than the Iceman that does. `I14677` and
`I14678` place at `G-PF3239` and `I15942` at `G-L91` on coverage, and **all three are excluded from
`G-Z6219` as well as `G-L166`**, which is stronger than the earlier prose, which only recorded the
`L166` exclusion. `UNTA58_68Sk1` returns no placement. Every rung of the ladder in
`iceman-y/README.md` now comes from a committed table through a committed script.

## 2026-07-27 — parked: open threads and two corrections

No analysis run. Session was review and course-correction; both corrections are recorded in place
(`NOTES.md`, `METHODS.md`) rather than only here.

**Corrections made**

- `NOTES.md`: the chrY-only emulation measured *recovery*, not *contamination*. `Z6215`/`Z6208`
  are not the same artifact as `FGC5687` — they show zero MQ0 and zero off-target under a
  whole-genome reference, so no read moved. The reverse simulation (tile reads from the paralogue,
  count how many land on chrY) has never been run and is the half that decides the question.
- `METHODS.md`: the "decoys can remove real signal" consequence is **aDNA-scoped**. The `-k101`
  collision filter plus a complete-Y reference plus a same-assembly collision index close the
  mechanism for modern saliva. Longer reads make the filter more conservative, not less.

**Open threads, roughly by value**

1. `iceman-y/prereg/Z6219_node.md` — registered outcome still reads "one site … no second
   population has been tested" and names P3 as unrun. Test C ran it at two further sites. The same
   amendment can discharge **F2**: across the five `Z6219`-derived Aesch/Muttenz men there is not
   one ancestral read anywhere on the backbone, and four are positively derived at `Z6043`
   (transversion). Caveat to state: `Z6043` is an `equivalent` marker on `G-PF3239`, so it
   establishes block membership rather than `PF3239` itself.
2. **Test B** — the co-segregation scan. Amendment 2 requires its marker set to come from YFull's
   32-SNP `L166` definition rather than `resources/marker_index.tsv.gz`; obtaining that list is a
   prerequisite and has not been done.
3. **Reverse mappability simulation** — see the NOTES entry above. Cheap, and it bounds a finding
   already committed.
4. **`Z6499` direction** — below vs parallel rests on parsimony alone; no sample derived at
   `Z6499` and ancestral at `L166` has been searched for.
5. **Protocol item D** (`Z6135`/`Z6209`) — post-hoc scan result, unregistered, needs its own
   document or an explicit decision to drop it.
6. **FC3 repair** — needs registration and a re-run; deliberately not patched into Test C.
7. **Furtwängler supplements via `openpyxl`** — check on the hand-rolled xlsx parser's numbers.
8. `annotate/y_novel_scan.py` — unbuilt. `WITH_HAPLOTYPERS=1` in `util/build_env.sh` untested past
   `bash -n`. `TODO.md` 1-5 unrelated mapping chores.
9. **Branch** — `tooling/se-adna-y-pipeline` sits well ahead of `master`, unpushed, with nothing
   else in flight. Recommended merging; not done.

**Self-control on a modern sample (planned, not started)**

The idea is a positive control: run the placement stack against a modern sample whose Y-SNP call is
already known from a commercial pipeline, at high coverage and with no deamination — the regime
where the call rules should be trivially right, which is what makes a failure informative. Design
notes that survive without reference to any particular sample:

- **What it does and does not exercise.** It tests `y_sites_pileup.py` -> `ylib.site_call()` ->
  `y_path_rank.py` end to end. It does **not** touch `ytree.py`: `markers/tree_local.tsv` is
  G-L166-only, so the tree machinery stays idle unless a second tree file is written for the
  sample's own branch. `y_tree_place.py --tree` already takes a path, so that costs nothing
  structurally; such a file names a living person's clade and belongs outside this repository.
- **Registrable prediction.** A commercial Y-only pipeline should disagree with a
  whole-genome-referenced call at X-homologous sites. The sites can be named in advance from
  `results/mappability/`, which makes this a real pre-registration rather than a comparison.
- **Reference mismatch is the practical obstacle.** Marker catalogues are GRCh38; a CHM13-based
  callset needs coordinate resolution first. No chain file is required — `y_mappability.py`
  discovers cross-build coordinates by modal implied position (`annotate/y_mappability.py:24`).
- **Archive the gVCF, not the VCF.** This method turns on ancestral and no-call states; a
  variants-only VCF discards exactly the rows that separate "tested and ancestral" from "never
  covered". `y_haplo_from_markers.py` reads gVCF for that reason.
- **Decoy cost measurement** (see METHODS.md) rides along on the same run at no extra cost.

Per-sample paths, identifiers and results stay outside this repository per AGENTS.md.

**Ideas parked**

- Custom metagenome decoy paper; eHOMD cherry-picking (`METHODS.md` records why the subset idea was
  rejected once already).
- Long-read data exists for at least one modern sample; **re-assembling parts of the Y scaffolding**
  from it is a separate project, noted so it is not lost.
- Broad/WARP vs nf-core/eager: still deferred, still awaiting a concrete trigger. `annotate/` is
  downstream of both and consumes BAMs and pileup tables, so nothing in flight is blocked.
