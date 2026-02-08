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
  | ./util/eager_repair_bam.py --strip-prefix --duplicate-merged --dup-orientation flip -o - --output-uncompressed \
  | samtools sort -n -@8 -m 2G -o stripped.qn.bam -

./util/eager_repair_bam.py --pair stripped.qn.bam -o paired.bam --output-uncompressed
```
Notes:
- `XQ:Z:M` marks original merged reads; `XQ:Z:D` marks synthetic mates.
- Pseudo-pairs enable insert-size modeling and PE aligner behavior, but do not create new sequence evidence.

## Remapping Workflow
Use the standard pipeline (e.g., `mapping/revert-bam.sh`) to create uBAMs and remap with BWA MEM.
For remap comparisons:
- Keep parameters constant across references.
- Compare mapping rates, mate discordance, and insert-size distributions.

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

## Metagenome Decoys
For Iceman, oral decoys produced sparse, spiky hits after human-reference filtering.
These reads are too few to affect human variant calling and are likely spillover from older reference mapping.
Recommendation: omit oral decoys for this dataset unless explicitly studying microbial signal.

## Interpretation Boundaries
Because the input is already human-filtered:
- Absence of strong microbial signal does not prove absence in the original sample.
- The practical conclusion is low impact on human analysis, not biological absence.
