# Methods (Ancient DNA Remapping Notes)

This document records practical methods and reasoning for re-mapping ancient DNA using the scripts in this repo. It is written for pre-filtered BAMs (e.g., reads already filtered to a human reference).

## Scope
- Applies to aDNA BAMs that have already been filtered against an older human reference.
- Examples reference the Otzi 2023 resequencing data (Cell Genomics 2023).
- Assumes reference builds produced by `mapping/GRCh38_bwa_index.sh`.

## Pair Reconstruction
When read pairs have been stripped, reconstruct pseudo-pairs for PE workflows:
```bash
samtools view -u input.bam \
  | ./util/eager_repair_bam.py --strip-prefix -o - --output-uncompressed \
  | samtools sort -n -@8 -m 2G -o stripped.qn.bam -

./util/eager_repair_bam.py --pair stripped.qn.bam --duplicate-merged --dup-orientation flip -o paired.bam --output-uncompressed
```
Notes:
- `XQ:Z:M` marks original merged reads; `XQ:Z:D` marks synthetic mates.
- Pseudo-pairs enable insert-size modeling and PE aligner behavior, but do not create new sequence evidence.
- `--pair` is the step that assigns READ1/READ2 flags for F/R. `samtools fixmate` does not create pairs.

## Remapping Workflow
Use the standard pipeline (e.g., `mapping/revert-bam.sh`) to create uBAMs and remap with BWA MEM.
For remap comparisons:
- Keep parameters constant across references.
- Compare mapping rates, mate discordance, and insert-size distributions.
Notes:
- `SANITIZE=true` in RevertSam(Spark) drops/normalizes problematic mates and can break F/R flag semantics. For EAGER-style data, prefer `SANITIZE=false` and skip MarkDuplicatesSpark.
- `samtools fastq` collapses multiple READ_OTHER records with the same QNAME down to one. If F/R reads lack READ1/READ2 flags, the R read disappears unless you set flags first.

## Insert-Size and QC
Use PE outputs for insert-size summaries:
```bash
samtools stats sample.bam | ./util/is_stats.py
```
A single dominant peak is preferred over multimodal peaks, especially when mixing libraries.

## Duplicate Handling
- MarkDuplicates on singletons over-marks duplicates; treat duplicate counts with caution.
- For pre-filtered aDNA, re-dedup can bias results when mates are missing.
- If needed, dedup only true PE reads or use DeDup on prefix-tagged data.
DeDup specifics:
- DeDup expects coordinate-sorted input and does not sort; output header is set to unsorted.
- DeDup removes duplicates (does not mark 0x400). Existing duplicate flags are preserved unless cleared first.
- DeDup requires `M_/F_/R_` prefixes unless `--merged` is used.

## DeDup Workflow (EAGER-style)
Recommended order for each library:
1) Restore prefixes and drop synthetic D mates; clear 0x400 if present.
2) Coordinate-sort.
3) DeDup on the library BAM.
4) Fix header or sort (DeDup writes `SO:unsorted`).
Then merge libraries and run a final DeDup on the merged BAM.
Notes:
- DeDup ignores read-group (RG). If library-specific dedup is required, split by RG first.
- Keep an “all contigs” coordinate-sorted BAM for later inspection; derive primary-only views from it.

## Fastq Splitting for PE/SE
If mates are missing from the BAM, do not rely on interleaved `bwa mem -p` with mixed PE/SE.
Preferred split:
- Use `samtools fastq -1/-2/-s` after READ1/READ2 flags are set.
- Route READ_OTHER to `-0` if you want to keep merged/SE reads.
Chunk boundary risk:
- `bwa mem -p` pairs only adjacent names and can split true pairs across chunk boundaries when singletons exist. Mitigations are separate PE/SE fastqs or larger `-K` to reduce splits.

## Metagenome Decoys
For Iceman, oral decoys produced sparse, spiky hits after human-reference filtering.
These reads are too few to affect human variant calling and are likely spillover from older reference mapping.
Recommendation: omit oral decoys for this dataset unless explicitly studying microbial signal.

## Interpretation Boundaries
Because the input is already human-filtered:
- Absence of strong microbial signal does not prove absence in the original sample.
- The practical conclusion is low impact on human analysis, not biological absence.

## Current Experiment Snapshot (Otzi D2049 Combined)
Dataset path used for baseline:
- `/mnt/AncientDNA/Iceman-2024/iceman.oetzi.UDG_D2049_combined.mapped_rmdup.bam`

Published BAM prefix composition:
- `M=138980507`
- `F=5545558`
- `R=3582337`

Reconstruction checkpoint (`pair.bam`) after `--strip-prefix -> qname-sort -> --pair --duplicate-merged --dup-orientation flip`:
- total records: `287088909`
- read1/read2: `144526065 / 142562844`
- singletons: `5182509`

Primary-assembly slice before DeDup (`pair.prim.bam`):
- total: `145064186`
- `M/F/R = 136687692 / 5138191 / 3238303`

Primary-assembly slice after DeDup (`pair.prim_rmdup.bam`):
- total: `138328921`
- `M/F/R = 131715629 / 3743710 / 2869582`

DeDup removals from log:
- merged removed: `4972063`
- forward removed: `1394481`
- reverse removed: `368721`
- total removed: `6735265`
- duplication rate reported by DeDup: `0.05` (rounded)
- exact from counts: `6735265 / 145035948 = 0.04644` (4.64%)

Observed loss rates by read class:
- `M`: `-4972063` (`-3.64%`)
- `F`: `-1394481` (`-27.14%`)
- `R`: `-368721` (`-11.39%`)

Interpretation used for this run:
- DeDup behavior is driven by `M/F/R + start/end/strand` rules and quality tie-breaking (not insert-size/template-aware duplicate logic).
- For this filtered aDNA dataset, keep DeDup output as one branch, but treat singleton-heavy unmerged duplicate calls as a workflow caveat.
