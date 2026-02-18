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

## Open Questions
- Quantify survivor vs removed read characteristics beyond prefix counts (MAPQ distribution, context near repetitive loci).
- Decide whether alt/decoy calls are retained only as technical appendix or excluded from variant interpretation entirely.
- Decide whether gVCF is required for all exploratory runs (skip when not needed).
