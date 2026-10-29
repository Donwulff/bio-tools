# Repository Situation Review

Reviewed **2026-10-28** against `master` at `43273c02ce4a0bd788118d4da6c2a46460894d29`,
the current scripts and committed public tables, the available session history, and GitHub PR
metadata. This is a documentation audit, not a new sequencing experiment or amendment to a
registered test. Modern sample identifiers, paths, genotypes and comparison counts are deliberately
excluded.

## Compact Handoff

Updated **2026-10-29** for the local checkpoint and reusable-tooling follow-up. The remote inventory below
retains its 2026-10-28 date; remote status was not rechecked as part of committing local work.

- The documentation reconciliation, reviewed Work research snapshots and marker-resolver
  corrections are preserved in local commits. No push or remote PR merge was performed.
- Subsequent comparator and shared-reference safeguards are preserved in focused local commits.
  These use 29 comparator, 27 reference and 41 resolver synthetic tests plus 40 tree checks;
  no sample analysis or registered outcome was rerun.
- Batch liftover now preflights every input, rejects unsupported gVCF confidence records from
  content and validates staged outputs before publishing a complete directory. This adds 26
  synthetic lifecycle/scanner tests; no production Picard conversion or sample analysis was run.
- Completed tooling: native gVCF block-base lookup, all-state quality thresholds, explicit filtered
  evidence and backward-compatible status-table consumption. Validation: 41 synthetic tests,
  40 tree checks and four public placement datasets reproduced without changed results.
- Recheck with `sh annotate/run_tests.sh` and `git diff --check`. Details are in
  [the resolver follow-up](RUNLOG.md#2026-10-28---native-gvcf-resolver-and-quality-filter-regressions)
  [the comparator/reference checkpoint](RUNLOG.md#2026-10-29---comparator-safeguards-and-shared-reference-validation)
  and [the batch-safety follow-up](RUNLOG.md#2026-10-29---ordinary-vcf-batch-liftover-safety).
- Comparator duplicate/allele safeguards, callable/expected-panel metrics and shared reference
  validation are implemented below. Next: resource/config/backend consolidation, panel provenance
  and AD/liftover-loss evidence. Full gVCF block liftover is separate
  and is **not** repaired by the native-query correction.
- Research decisions remain in the open-work table below. Do not relaunch experiments, change
  registrations or use private sample data as an implicit continuation of a tooling fix.
- The commit clarification applied to the previous checkpoint; the subsequent request resumed
  reusable-tooling work. Research extensions still need their own decisions and registrations.
- Standing user instruction: commit finished, validated batches automatically; do not push or
  merge remote PRs unless requested. This checkpoint does not implement the remaining research tasks.

## Start Here

| Question | Record to use |
|---|---|
| What is the current Iceman Y interpretation? | [Case verdict](iceman-y/README.md), with the qualifications below |
| Which observations and placements are reproducible? | [Committed public results](iceman-y/results/), [local tree](markers/tree_local.tsv), [run ledger](RUNLOG.md) |
| How were the selected EAGER BAMs reconstructed and remapped? | [Methods](METHODS.md), [pairing and DeDup troubleshooting](NOTES.md) |
| What goes wrong with hs1/GRCh38 liftover and comparisons? | [Liftover guide](annotate/LIFTOVER.md) |
| Which scripts implement each operation? | [Annotation tools](annotate/README.md), [mapping tools](mapping/README.md), [utilities](util/README.md) |
| Which tests were registered, and which remain open? | [Registrations](iceman-y/prereg/), [extension protocol](iceman-y/PROTOCOL_extending_analyses.md), open-work table below |
| Where are the newer Work assessments? | PR inventory below and the locally preserved [candidate assessment](iceman-y/research/otzi_literature_candidates.md) / [priorities](iceman-y/research/priorities_and_source_watch.md) |
| Can the workflow be kept in a notebook? | [Jupyter template](notebooks/otzi_reanalysis_template.ipynb); templates are not evidence of execution |

`NOTES.md` and `RUNLOG.md` are chronological working records. An early statement is not current
merely because it appears first. Later corrections, dated outcomes and their tables take precedence;
this review points to those corrections without rewriting preregistrations after the fact.

## Branches and ChatGPT Work

GitHub was checked on the review date. The working tree was clean at entry. Local
`tooling/se-adna-y-pipeline`, `master`, the locally recorded `origin/master` and GitHub's live
`master` all point to `43273c0`; the July note saying that tooling branch was ahead and unmerged is
now historical.
No merge, rebase, commit or push was performed in this review.

| Work item | Verified state | Content and consequence |
|---|---|---|
| [PR #7](https://github.com/Donwulff/bio-tools/pull/7), `research/otzi-literature-candidates-2026-10-27` | Open; head `2663d4c74cb8ce1946d909caa65c55e2ddd2f1f3` | Literature assessment and accession TSV; no reads downloaded or genotypes tested in that work |
| [PR #8](https://github.com/Donwulff/bio-tools/pull/8), `docs/otzi-research-value-and-source-watch` | Open; head `fd7264d69dc72d73f61ec7e93ec95536dd1df323` | Research-value assessment and source-watch description; no changed placements or registered outcomes |

The two PRs add different documents from the July baseline; neither includes the other's full
assessment. Their exact reviewed source versions are
[candidate research](https://github.com/Donwulff/bio-tools/blob/2663d4c74cb8ce1946d909caa65c55e2ddd2f1f3/iceman-y/research/otzi_literature_candidates.md),
[accession crosswalk](https://github.com/Donwulff/bio-tools/blob/2663d4c74cb8ce1946d909caa65c55e2ddd2f1f3/iceman-y/research/otzi_literature_candidates.tsv), and
[priorities/source watch](https://github.com/Donwulff/bio-tools/blob/fd7264d69dc72d73f61ec7e93ec95536dd1df323/iceman-y/research/priorities_and_source_watch.md).
The continuation of this review copied all three documents unchanged into `iceman-y/research/`:
[candidate assessment](iceman-y/research/otzi_literature_candidates.md),
[accession crosswalk](iceman-y/research/otzi_literature_candidates.tsv), and
[priorities/source watch](iceman-y/research/priorities_and_source_watch.md). Source commits are
the PR heads above. This is local documentation integration, not a GitHub PR merge or a new
experimental result. Their source-date claims and inherited terminology are retained as historical
assessment text; the qualifications in this review and the current verdict govern interpretation.

PR #7's bounded shortlist is XAN027 (Y-capture manifest and filtering provenance first),
SEC002/I18968 (new versus previously merged library identity and target overlap), SP193 (shotgun
coverage gate), and I13517 (damage-selected data and locus callability). It also traces additional
Iceman tissue/water metagenomic deposits under PRJEB94382; retention of authentic host-Y reads is
unverified. None of these leads establishes a new branch, independent patriline or usable marker
depth. The PR records access gaps and stopping conditions rather than treating unavailable data
as absent. Its maternal literature update is a separate question from paternal placement.

PR #8 usefully separates value for Iceman, value for reusable tooling and feasibility. It reports a
monthly external source watch. The imported document preserves that report; this checkout contains
no verified scheduler configuration or watch
execution log; its current operation cannot be inferred from that prose. No separate completed
Z6499 discriminator-search result was located among the available PRs or committed tables.

## What Is Already Preserved

### Remapping and DeDup

The public D2049 checkpoints and final four-library totals are already recorded in
[METHODS.md](METHODS.md#current-experiment-snapshot-otzi-d2049-combined) and
[RUNLOG.md](RUNLOG.md#confirmed-counts):

| Checkpoint | Records |
|---|---:|
| Original D2049 BAM | 148,108,402 |
| Original merged `M_` records | 138,980,507 |
| Reconstructed input, including one synthetic mate per merged record | 287,088,909 |
| Primary alignments after remapping | 287,088,909 |
| D2049 primary-assembly, primary-alignment input to DeDup | 145,064,186 |
| D2049 output from DeDup | 138,328,921 |
| D2049 records removed by DeDup | 6,735,265 |
| Four-library merge before final DeDup | 668,715,804 |
| After final DeDup | 656,719,414 |
| Removed in final merged pass | 11,996,390 |

The reconstruction identity is `148108402 + 138980507 = 287088909`. Synthetic mates must be
removed before comparing with original molecule/read counts. The D2049 DeDup removals break down
as M=4,972,063, F=1,394,481, R=368,721; the name-join comparison reproduced the log. The reduction
between the original BAM and the final per-library output also includes assembly/alignment
selection, so that entire difference cannot be labelled duplicates.

Lessons retained from this session: establish F/R read-number flags before `samtools fastq`, disable
RevertSam sanitizing when preserving missing-mate records is required, restore prefixes/drop D
mates before DeDup, distinguish name grouping from coordinate order, and coordinate-sort actual
DeDup output rather than simply claiming it sorted in its header. Per-library deduplication followed
by merged deduplication was intentional; the algorithm's singleton behavior remains a caveat.

### DeepVariant

Both small-model-on and small-model-off Iceman calls completed. The committed
[A/B summary](iceman-y/results/iceman_smallmodel_compare.summary.tsv) gives 38,024 shared-record
GT differences out of 5,010,771 (0.758845%) over all records and 671 out of 3,969,848 (0.016902%)
over shared PASS records. M201, L91, L166 and Z6208 retained their genotypes in the
[marker comparison](iceman-y/results/iceman_smallmodel_compare.y_markers.tsv).

These are concordance measurements, not measured error rates or proof that either model is better
on aDNA. The actual runs differed in region scope; the on run traversed zero-read non-primary
contigs for gVCF output. The correction in
[RUNLOG.md](RUNLOG.md#deepvariant-ab-plan-small_model-onoff) distinguishes executed commands from
the cleaner templates that were initially documented. No matched modern truth-set benchmark was
completed in the available record.

The long CPU-only inference and protobuf-heavy postprocessing did finish. Persisting intermediates
allows manual restart at stage boundaries; the interrupted run was not demonstrated to resume
inference midway through a partially written TFRecord. Container storage, logs and process state
must be inspected before starting a duplicate run. Output compression must be checked from content,
since the original bare `.vcf` and `.gvcf` names contained compressed data.

### Y Interpretation and QC

February's label ranking established derived markers but did not establish their phylogenetic order.
The July read-level work supports the project's **G-L166 placement**. Z6208 remains derived; its
treatment as a distinct downstream terminal node is rejected in the **local tree**, not its allele
call. The provisional PF3239 -> Z6219 -> L166 structure and the proposed Z6499 split are recorded
with statuses and limitations in the [case verdict](iceman-y/README.md).

The committed local tree covers this G subtree, not every Y haplogroup. A modern control on another
branch needs its own verified topology kept with the external analysis. February's preliminary
VCF comparison does not complete the July proposal to validate the newer read-level placement stack.

The seven positively placed Swiss men occur at three sites; that is not seven demonstrated
independent paternal lineages. I5118 has sparse supporting transversions, one molecule per tested
site. The placement tables also return Aesch7 at L166 with `weak_ancestral=1;single_marker`, so a
claim that I5118 is literally the only other tool-returned L166 placement needs that qualification.
Z6499's below-versus-parallel direction remains a parsimony hypothesis. An unresolved child or
no-call cannot establish a universally basal terminal-star assignment.

An expected derived marker called confidently ancestral is a useful discordance to investigate;
no-call measures missing evidence, which is especially common in pre-filtered aDNA. Neither is
automatically a caller error: recurrence, catalogue mistakes, contamination, paralogues, low depth
and reference/filter bias must remain distinguishable. `DP` describes depth, `AD` allele depth;
the ranker's DP sums do not count independent supporting molecules. Equivalent SNP aliases and
overlapping result files can otherwise inflate apparent support.

### Liftover and Modern Controls

The session includes completed modern VCF/gVCF conversion attempts and a preliminary comparison
with direct GRCh38 calls. Per-sample evidence belongs outside this repository. What is preserved
here is the method and its limitations in [annotate/LIFTOVER.md](annotate/LIFTOVER.md).

The earlier recommendation to trust direct GRCh38 automatically was not established by that
comparison. Three callsets made from the same reads are correlated. Concordance at familiar Y
markers is a useful control, while agreeing no-calls over millions of catalogue records is not a
general accuracy score. The VCF and gVCF marker panels and callable denominators must be matched.
Picard's treatment of gVCF blocks is an implementation limitation, not just a difficult reference
transition. Consolidation remains partial: validation/naming now share one implementation across
marker queries and shell preparation; reference discovery and backend execution remain duplicated.

The subsequent 2026-10-28 code follow-up adds indexed native-reference lookup and all-state
quality filtering to the marker resolver, with synthetic regression tests. Drivers use that same
engine; old sample outputs were not regenerated. This does not repair Picard gVCF block conversion,
validate the preliminary modern comparison, or alter registered read-level outcomes.

The 2026-10-29 follow-up adds comparator duplicate/allele guards, callable-only agreement, explicit
expected-panel QC and safe reference staging/validation. Default reference checks are metadata and
small probes; full M5/byte-layout verification is opt-in. Synthetic tests validate these operations,
not any sample placement or the correctness of a supplied expected panel.

## Open Work, Without Reopening Completed Experiments

### Registration Reconciliation

This is a dated reading of existing records, not an amendment to their predictions or outcomes.
The original [Z6219 registration](iceman-y/prereg/Z6219_node.md) remains intact alongside the
separate [Test C registration and outcome](iceman-y/prereg/testC_aesch_muttenz.md).

| Question | What the existing records establish | What remains unresolved |
|---|---|---|
| Original P3 / Test C PC2: a split at Z6219 | Test C reports all five callable individuals derived, with no split; uniformity was explicitly registered as uninformative | Do not mark the original split prediction fulfilled by Test C's PC1 result |
| Test C PC1: intermediate state beyond Oberbipp | The recorded outcome reports four individuals callable at both markers, three documented kin groups at two additional sites | This is site replication, not demonstrated independent patrilines; FC3 remains not evaluable |
| Original F2: Z6219-derived but PF3239-ancestral or outside G | Test C reports no ancestral backbone calls and derived Z6043 in four of the five Z6219-callable men; their direct PF3239 rows are two uncovered and three one-read damage-prone no-calls | Z6043 is an equivalent block marker, not a direct PF3239 call; these no-calls cannot be turned into a completed per-person F2 test |
| Test B source membership | The committed marker list contains 22 resolved positions and the existing scan used those positions | Its title refers to 32 source groups; the complete authoritative source roster and ten unresolved groups are not reconciled by that title or a later live tree |

Any new independent-sample study, FC3 repair or protocol-D validation needs its own prospective
decision rules. A documentary cross-reference cannot supply the missing observation or repair a
rule after its outcome is known.

### Remaining Decisions

| Item | Current bounded status / next action |
|---|---|
| Work PRs #7 and #8 | Source documents preserved locally at the reviewed commits; remote merge/disposition remains pending, with no changed results |
| Z6219 registration P3/F2 | Documentary reconciliation above preserves both outcomes; any formal dated addendum must retain their distinct predictions and missingness |
| Test B source membership | 22-position execution exists; reconcile the authoritative 32-group source roster and ten unresolved groups before claiming an exhaustive test |
| Reverse mappability | Unrun: simulate non-Y/paralogue-origin reads accepted on Y; locus-origin recovery alone does not measure this |
| Z6499 direction | No completed search for derived Z6499 / ancestral L166 found; inspect the separate Work result if one becomes available before duplicating it |
| Z6135/Z6209 (protocol D) | Post-hoc observation; prospective validation or explicit disposition remains required |
| FC3, supplement parsing | FC3 repair needs a separate registration/rerun; independent spreadsheet-parser check remains unrecorded |
| Modern three-way QC | Repeat with identical marker definitions, native gVCF evidence, explicit missingness and trusted expected-path definitions; keep every result external |
| Liftover tooling | Reference checks and ordinary-VCF batch preflight/staged publication implemented; resource/config/backend consolidation and other drivers' publication safety remain open, as does gVCF block conversion |
| Marker comparator | Duplicate/allele guards and callable/expected-panel reports implemented; verified crosswalks, AD and exact reject provenance remain open |
| Marker resolution/quality | Native-reference block lookup and all-state thresholds implemented with synthetic tests; matched-panel sample QC, panel hashes and broader allele validation remain pending |
| Ranking output | Separate marker evidence from clade placement; retain current heuristic weights as descriptive, not learned probabilities |
| Novel scan and external callers | `y_novel_scan.py` is still absent; PathPhynder/Yleaf cross-checks and optional installs have no recorded end-to-end validation |
| Broader analyses | Exomiser/PharmCAT, phenotype review and multi-platform modern integration were discussed, not completed in this record |

The passed tree tests check code and the local topology's encoded rules, not the truth of every
scientific assumption. This review ran the existing 40 checks and reproduced placements from the
committed public pileups in temporary storage; it did not remap samples, change thresholds or
rerun claim-bearing experiments. See the dated [review ledger](RUNLOG.md#2026-10-28--situation-review).
