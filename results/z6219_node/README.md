# Outputs of Tests A and B, `PREREG_Z6219_node.md`

Run 2026-07-26. Commands in `RUNLOG.md`, findings in `NOTES.md`.

## Which directory is canonical

`results/swiss15/` and `results/unhedged/` remain the canonical genotype tables. They are not
superseded by anything here.

`results/testB/` and `results/testB_unhedged/` are the full output of one `y_genotype_batch.sh`
invocation each. They exist because the YFull 32-SNP marker set needed genotyping and the batch
script processes every marker set in `markers/` in one pass — so re-running it for Test B
regenerated the older tables as a by-product. Those by-products were kept rather than deleted
because the `RUNLOG.md` command reproduces the whole directory, and a documented command that
produces files not present in the repository is a command nobody can check.

All ten regenerated tables were verified **byte-identical** to their canonical counterparts. That
is the point of keeping them: they are a reproducibility regression with the diff already run, not
a second opinion. If they ever stop matching, something changed that should not have.

The genuinely new tables are `*_yfull_L166_defining.tsv` (22 testable positions of YFull's 32-SNP
`L166` definition) and `swiss_yfull_pooled.tsv` (Family A as one observation).

## Files here

| file | what it is |
|---|---|
| `cosegregation_yfull32.tsv` | Test B output: per-marker pattern across upstream cohort vs outgroup |
| `I5118_read_evidence.txt` | Test A: per-read evidence at `Z6219` and `L166`, plus the library damage profile |
| `familyA_backbone_pooled.tsv` | Family A pooled backbone calls, for the F2 haplogroup-G check |
| `iceman_yfull_L166_defining.tsv` | E1: Iceman genotyped at all 22 testable YFull `L166`-defining positions |
| `iceman_E1_classification.tsv` | E1: the same table under the counting rule registered before the run |
| `iceman_E1_read_evidence.txt` | E1: library damage profile, plus per-read evidence at `FGC5671` and `Z6499` |
| `e1_mappability_{anc,der}.tsv` | E1: allele-aware recovery at `FGC5671` and `Z6499`, 35/45/60/100 bp |

## Reading `iceman_E1_classification.tsv`

The 9 `control(ok)` rows are `markers/L166_defining.txt`, already known derived in the Iceman before
E1 was registered. They are excluded from the counts and exist only to invalidate the run if they
are not 9/9.

`untested` again means **not tested**. All five untested positions fail the inherited 30% MQ0
rejection threshold or have no coverage; three of them carry `DERIVED`-direction calls that the
threshold discards. Those three are not evidence of anything and must not be promoted by relaxing
the rule after the fact.

## Reading `cosegregation_yfull32.tsv`

`uninformative` means **not tested** — no coverage, a single damage-prone read, or mixed. It is not
a negative result and must not be counted as one. F3 in the pre-registration requires markers tested
and found not to co-segregate; 17 of the 22 positions never reached that bar.
