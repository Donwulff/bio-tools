# Notes

**Dataset**
Cell Genomics 2023 short article (Sep 13, 2023): "High-coverage genome of the Tyrolean Iceman reveals unusually high Anatolian farmer ancestry." This note set refers to the Otzi 2023 resequencing BAMs.

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
