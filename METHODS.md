# Methods (Ancient DNA Remapping Notes)

This document records practical methods and reasoning for re-mapping ancient DNA using the scripts in this repo. It is written for pre-filtered BAMs (e.g., reads already filtered to a human reference).

## Scope
- Applies to aDNA BAMs that have already been filtered against an older human reference.
- Examples reference the Otzi 2023 resequencing data (Cell Genomics 2023).
- Assumes reference builds produced by `mapping/GRCh38_bwa_index.sh`.

## Data Governance
- Published ancient/public data can be documented and versioned in this repository.
- Modern/private sample analysis must be performed in external paths and must not be committed with sample identifiers or per-sample variant outputs.
- Keep committed content to reusable scripts, parameter templates, and de-identified method notes.

## Path Conventions
- `<raw_data_dir>`: location of original input BAM files.
- `<analysis_dir>`: location of per-run working/output files.
- `<legacy_data_dir>`: location of legacy comparison VCF datasets.

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

### Supporting references (oral/eHOMD decoy rationale)
The decoy practice itself is inherited from the general hs38d1/GRCh38DH convention; no
publication was found that constructs and evaluates an *oral-specific* decoy set. The
supporting literature establishes the problem and recommends the remedy rather than testing it:

- Samson, C. A., Whitford, W., Snell, R. G., Jacobsen, J. C. & Lehnert, K. Contaminating DNA
  in human saliva alters the detection of variants from whole genome sequencing.
  Sci. Rep. 10, 19255 (2020). — establishes that non-human saliva DNA aligns to the human
  reference under standard methodology and alters resulting genotypes.
- Kumar, A. et al. Microbial contamination and composition of oral samples subjected to
  clinical whole genome sequencing. Front. Genet. 14, 1081424 (2023). — characterises the
  contamination and recommends adding prevalent oral microbial genomes as decoy. Note it
  **recommends only**; its own analysis used plain GRCh38.p13 with non-human reads filtered out.
- Chrisman, B. et al. The human "contaminome": bacterial, viral, and computational
  contamination in whole genome sequences from 1000 families. Sci. Rep. (2022).
  doi:10.1038/s41598-022-13269-z, PMC9198055. — **the important caveat for this repo.**
  Y-chromosome fragments absent from GRCh38 commonly mismap to bacterial reference genomes
  (only about half of chrY is in GRCh38, and it is repeat-rich), producing >100 bacterial
  species falsely associated with sex. Constructs no decoy; recommends masking the
  problematic k-mers instead.

Consequence for Y work: with bacterial decoys present, genuine chrY reads from poorly
catalogued regions have somewhere wrong to go, so decoys can *remove real signal*, not just
noise. This is an independent, published account of the same artifact class documented in
NOTES.md (MAPQ-0 pileups, the chrY:11.1-11.7 Mb no-go zone, collapsed repeats near 56.7 Mb),
and it is a second reason to keep oral decoys off for the Iceman analysis.

### Read pre-filtering vs decoy-in-reference
An, Z., Cha, J. H., Lee, K. H. & Lee, I. Metagenome-assembled genomes enhance bacterial read
decontamination and variant calling in oral samples. iScience (14 Oct 2025), PMC12616088 —
reports that MAG-augmented catalogues (HROM, 72,641 genomes) beat eHOMD on 5 of 6 SNP metrics.
Two caveats before reading that as a verdict on this pipeline:

- **Different mechanism.** They use Kraken2 pre-alignment classification and discard reads
  classified as bacterial; they are not adding decoys to the reference. They explicitly
  position this against the decoy approach. Their eHOMD arm still improved over no
  decontamination at all — the comparison is better-vs-good, not good-vs-harmful.
- **Pre-filtering is the riskier design.** Discarding reads before alignment makes the
  classifier the final arbiter, with no scoring contest and no surviving read to audit.
  Where a bacterial genome is homologous to a human region, dropout can be allele-specific
  and manufactures reference bias invisibly. Their filtering was for assembly quality
  (completeness, contamination, chimerism), which is orthogonal to human-homology collisions.
  The `-k101` collision filter in `mapping/GRCh38_bwa_index.sh` targets exactly that and has
  no equivalent in their pipeline.

Version note: their eHOMD arm filtered to 8,067 genomes, consistent with the V10.1+ generation
(8,622 genomes, 2023-03-18) rather than the ~2,123-genome pre-10.1 releases. `VERSION_ORAL`
defaults to 9.15, i.e. the smaller older generation — which is also why 10.01+ needs 64 GB.

### Idea not pursued: eHOMD subset by hit count
Selecting only the eHOMD genomes that actually attract reads would cut decoy size and memory
substantially. Rejected as not worth it: determining which organisms are *universally* common
would require mapping or BLASTing a large, internationally representative set of saliva-derived
sequences, and the resulting subset would be a judgement call that needs revisiting per
population and per eHOMD release. Using the whole collection with human-analogous sequences
removed is simpler and has no selection bias to defend.

## Interpretation Boundaries
Because the input is already human-filtered:
- Absence of strong microbial signal does not prove absence in the original sample.
- The practical conclusion is low impact on human analysis, not biological absence.

## Current Experiment Snapshot (Otzi D2049 Combined)
Dataset path used for baseline:
- `<raw_data_dir>/iceman.oetzi.UDG_D2049_combined.mapped_rmdup.bam`

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

## Final Merge Snapshot (D2049-D2052)
Merged primary-assembly BAM before final DeDup (`...sort.bam`):
- total: `668715804`
- mapped: `668606112` (99.98%)
- paired in sequencing: `649228301`
- read1/read2: `639093736 / 10134565`
- properly paired: `641169636` (98.76%)

After final merged DeDup (`...sort_rmdup.bam`):
- total: `656719414`
- mapped: `656609722` (99.98%)
- paired in sequencing: `639445514`
- read1/read2: `630194482 / 9251032`
- properly paired: `631748233` (98.80%)

Observed final-pass delta:
- removed reads: `11996390` (`1.79%`)
- mapping rate remained stable (`99.98%`)

## Removed-Read Characterization (First Library)
For D2049, an exact removed SAM/BAM was created by join-diff between name-sorted pre/post DeDup BAMs.
Important implementation details:
- use tab-delimited join: `join -t$'\\t'`
- use name order compatible with join (`samtools sort -N` or explicit `LC_ALL=C sort`)

Reference command pattern:
```bash
samtools view -H pre.Nsort.bam > removed.sam
join -t$'\\t' -v1 <(samtools view pre.Nsort.bam) <(samtools view post.Nsort.bam) >> removed.sam
samtools view -b -o removed.bam removed.sam
samtools flagstat removed.bam
```

Observed removed BAM (`...pair.prim_rmdup.Nsort.removed.bam`) summary:
- total removed records: `6735265`
- paired in sequencing: `5768582`
- read1/read2: `5492819 / 275763`
- properly paired: `5658405` (98.09%)

Interpretation:
- removals are concentrated in read1-heavy classes, consistent with DeDup rule-based treatment of `M/F/R` in singleton-heavy data.

## DeepVariant Interrupt/Resume
DeepVariant `run_deepvariant` is not a checkpoint manager by itself. Resume depends on intermediate files surviving on host storage.

Recommended run option:
```bash
--intermediate_results_dir=/output/dv_tmp_<run_id>
```

Resume rules:
- If `call_variants_output*.tfrecord.gz` is complete: run `postprocess_variants` only.
- If only `make_examples.tfrecord@N.gz` exists: run `call_variants`, then `postprocess_variants`.
- If intermediates are missing (container `/tmp` lost): rerun full workflow.

Manual resume pattern:
```bash
# 1) call_variants from saved make_examples outputs
/opt/deepvariant/bin/call_variants \
  --examples /work/dv_tmp/make_examples.tfrecord@8.gz \
  --outfile /work/dv_tmp/call_variants_output.tfrecord.gz \
  --checkpoint /opt/models/wgs

# 2) postprocess_variants
/opt/deepvariant/bin/postprocess_variants \
  --ref /work/ref.fa \
  --infile /work/dv_tmp/call_variants_output.tfrecord.gz \
  --outfile /work/sample.vcf \
  --gvcf_outfile /work/sample.gvcf \
  --nonvariant_site_tfrecord_path /work/dv_tmp/gvcf.tfrecord@8.gz \
  --checkpoint_json /opt/models/wgs/model.example_info.json \
  --cpus 8
```

If original run used `--disable_small_model=false`, include:
```bash
--small_model_cvo_records /work/dv_tmp/make_examples_call_variant_outputs.tfrecord@8.gz
```

## Single-End aDNA from Raw FASTQ (no pre-filtered BAM)

The pair-reconstruction workflow above exists to undo read merging in already-mapped,
pre-filtered BAMs. It does not apply to raw single-end FASTQ obtained from SRA/ENA: there are
no pairs to reconstruct and no `M_/F_/R_` prefixes to restore. `revert-bam.sh` and
`eager_prep.sh` are skipped entirely for this input class.

Entry points:
- `annotate/fetch_ena_runs.sh <accession> <outdir> [alias_regex]` — writes a manifest and
  stages runs with MD5 verification against ENA-published checksums; re-runnable.
- `mapping/map_se_adna.sh <fastq.gz> <sample_id> [run_id] [threads]` — `bwa aln`/`samse`,
  coordinate sort, `samtools markdup -r`, and a chrY coverage summary.

Aligner choice: `bwa aln`/`samse`, **not** `bwa mem`. `bwa mem` is designed for reads >=70 bp
and loses sensitivity below that; aDNA capture reads are typically 35-90 bp. Parameters follow
nf-core/eager 2.5.0 defaults (`-n 0.01 -k 2 -l 1024`); `-l 1024` exceeds any realistic aDNA read
length and therefore disables seeding, so the whole read including its damage-bearing termini
participates in alignment.

Two consequences to record whenever this path is used:
- `bwa aln` is **not ALT-aware**. `bwa-postalt.js` is not applied and the reference's `.alt`
  companion goes unused. Acceptable for Y marker genotyping, where MAPQ filtering does the
  equivalent work, but it is a real difference from any `bwa mem` run on the same reference.
- Deduplication is **mandatory**, not optional. Capture libraries are PCR-amplified, and any
  rule that counts "independent reads" is void if duplicate copies of one molecule inflate
  apparent support.

`samtools markdup` is usually documented with a `collate` / `fixmate -m` / `sort` prerequisite,
which is a **paired-end** requirement: `fixmate -m` exists to add mate-score tags. On single-end
input `markdup` works directly off the coordinate-sorted BAM, and `map_se_adna.sh` calls it that
way. Verified rather than assumed — six synthetic reads (three forward at one position, two
reverse at the same position, one elsewhere) give `SINGLE: 6`, `DUPLICATE SINGLE: 3`,
`WRITTEN: 3`, keeping one read per strand-position group. Note it is strand-aware: the forward
and reverse reads sharing a start coordinate are correctly *not* duplicates of each other.

It does emit `warning, unable to calculate estimated library size ... Read pairs 0`. That is
cosmetic on this input class — `PAIRED: 0` is the expected state for single-end data, and the
library-size estimator is a paired-end statistic. Do not read it as a dedup failure.

Do not pursue identical alignment steps across incompatible library types for their own sake.
Comparability between a paired-end shotgun dataset and a single-end capture dataset lives
downstream of alignment — same reference, same marker set, same pileup thresholds — not in a
shared script name.

## Y Haplogroup Marker Workflow (Current)
For GRCh38-based Y haplogroup checks:
1. Merge sample chrY calls with YBrowse hg38 marker VCF (`snps_hg38.vcf.gz`).
2. Extract derived calls against marker ancestral allele (`INFO/AA`), with caller-specific site filtering.
3. Score candidate clades across datasets via support/conflict counts (not only "deepest label string").

Implementation in repo:
- `annotate/y_haplo_from_vcf.sh`
- `annotate/y_haplo_from_markers.py`
- `annotate/y_clade_consistency.py`
- `annotate/y_path_rank.py`
- `annotate/y_markers_pileup.py` — pileup of *named* catalogue markers, resolved via
  `resources/marker_index.tsv.gz`.
- `annotate/y_sites_pileup.py` — pileup of bare coordinates, for sites that have no catalogue
  name. Uncatalogued ("novel") candidate sites are exactly this case, so they cannot be queried
  by the marker script. Emits `pct_mq0` and `n_mq60` alongside allele counts.

### MAPQ auditing is not optional
Depth alone cannot distinguish a real variant site from a collapsed-repeat pileup. Sites with
entirely unremarkable depth (DP 4-12) have been rejected at 39-88% MAPQ-0 reads, including
transversions that had previously been presented as the strongest damage-immune evidence
available. Any claim about a chrY site must report the MAPQ distribution at that site, which is
why `y_sites_pileup.py` emits those columns unconditionally rather than on request.

DeepVariant-specific setting:
- use `--site-filter-mode deepvariant` in `y_haplo_from_vcf.sh` (`FILTER=PASS || RefCall`).

Attribution:
- Root-to-tip style support ranking and damage-aware SNP accounting follow published ancient-DNA practice:
  - https://doi.org/10.1101/2024.03.13.584607
- Parsimonious ancestral/derived path interpretation is in the same spirit as PathPhynder usage in that workflow.

### Optional external haplotypers

The placement logic here is an independent implementation of the approach cited above; PathPhynder
itself has never been run against these samples. Two external callers are available as cross-checks
and are installed by `util/build_env.sh` only when `WITH_HAPLOTYPERS=1`:

- **pathPhynder** (`ruidlpm/pathPhynder`, R) with **phynder** (`richarddurbin/phynder`, C). An
  independent placement of the same BAMs, against its own reference tree. Its R dependencies
  (`optparse`, `phytools`, `scales`, `ape`) are all packaged on Ubuntu; only phynder needs building,
  against the static htslib the same script produces.
- **Yleaf** (`genid/Yleaf`, Python, source only -- not on PyPI). Relevant because the
  `haplogroup.info` "All Ancient DNA" compilation states Yleaf v2 as its caller, so this is the tool
  whose output the disputed `G-L166*` labels derive from. **Use the v2 tag when comparing against
  those labels.** The current release is v4.x, and on v4 any disagreement is equally well explained
  by four major versions of tree updates -- which would defeat the purpose, since the question is
  whether the automated caller or the hobbyist re-derivation layered on top introduced the error.

Neither is invoked by any script under `annotate/` or `mapping/`. A tool earns unconditional install
logic when a committed pipeline script calls it; until then it is recorded, not installed.

Current result snapshot (Iceman):
- Illumina branch supports `M201`, `L91`, `L166`, and `Z6208` as derived.
- SOLiD branch remains conflict-heavy/inconclusive for terminal placement with current panel.

## Dual Liftover Workflow Checklist (Modern Y, hs1 <-> hg38)
Use this when validating both directions:
1. `hg38 markers -> hs1` (score on native hs1 sample).
2. `hs1 sample -> hg38` (score on hg38 markers/resources).

Practical guardrails:
- Prefer non-extended references for this task:
  - hs1: `chm13v2.0_maskedY_rCRS.fa(.gz)`
  - hg38: `GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna`
- Ensure reference sidecars exist before liftover:
  - `<ref>.fai`
  - dictionary (`.dict`), accounting for compressed FASTA naming (`ref.dict` vs `ref.fa.dict` aliasing).
- Sanitize marker VCF before liftover (`ALT!=REF`, no duplicate ALT alleles).
- Keep private modern outputs outside repo; commit only scripts/method text.
