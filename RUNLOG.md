# Run Log (Otzi Reanalysis)

This file is the operator-facing ledger for reproducibility. Keep methods in `METHODS.md`, interpretation in `NOTES.md`, and put exact run invocations/results here.

## Scope
- Dataset: Otzi 2023 resequencing BAMs, EAGER-style prefixed reads (`M_/F_/R_`), published in filtered form.
- Goal: reconstruct PE semantics where possible, remap, evaluate duplicate handling, and compare variant-calling branches.

## Canonical Inputs
- Baseline BAM: `/mnt/AncientDNA/Iceman-2024/iceman.oetzi.UDG_D2049_combined.mapped_rmdup.bam`
- Reference (primary run): `mapping/index/hg38p14DH3630O.fa`

## Key Decisions Locked
- Keep `mapping/revert-bam.sh` stage-gated (`ENABLE_ALIGN`, `ENABLE_MARKDUP`).
- Use explicit `samtools fastq` split stream before `bwa mem -p` for mixed PE/SE handling.
- For DeDup, use primary-assembly slice for production branch.
- Treat DeDup output as one branch; singleton-heavy over-pruning is a documented caveat.

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
bcftools stats /home/jsantala/iceman.vcf > /home/jsantala/iceman.vcf.stats
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
  -v /mnt/mirrored/iceman-reanalysis/dedup_out50:/input \
  -v /home/jsantala:/output \
  google/deepvariant:1.10.0 /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/output/src/bio-tools/mapping/index/hg38p14DH3630O.fa \
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
  -v /mnt/mirrored/iceman-reanalysis/dedup_out50:/input \
  -v /home/jsantala:/output \
  google/deepvariant:1.10.0 /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/output/src/bio-tools/mapping/index/hg38p14DH3630O.fa \
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
bcftools stats /home/jsantala/iceman_sm_on.vcf.gz  > /home/jsantala/iceman_sm_on.vcf.stats
bcftools stats /home/jsantala/iceman_sm_off.vcf.gz > /home/jsantala/iceman_sm_off.vcf.stats

bcftools isec -p /home/jsantala/iceman_sm_cmp \
  /home/jsantala/iceman_sm_on.vcf.gz /home/jsantala/iceman_sm_off.vcf.gz

wc -l /home/jsantala/iceman_sm_cmp/0000.vcf \
      /home/jsantala/iceman_sm_cmp/0001.vcf \
      /home/jsantala/iceman_sm_cmp/0002.vcf \
      /home/jsantala/iceman_sm_cmp/0003.vcf
```

## Open Questions
- Quantify survivor vs removed read characteristics beyond prefix counts (MAPQ distribution, context near repetitive loci).
- Decide whether alt/decoy calls are retained only as technical appendix or excluded from variant interpretation entirely.
- Decide whether gVCF is required for all exploratory runs (skip when not needed).
