# Pre-registration: L166 status of Swiss Late Neolithic G2a2a individuals

Written 2026-07-25, **before any read from this dataset has been examined**. Hypotheses,
decision rules and declared limitations below are fixed as of that date. Results go in
`RUNLOG.md` and `NOTES.md`; this file is not to be edited to match them.

## Background

The Tyrolean Iceman is terminal at `G-L166` by read-level evidence
(`results/iceman_y_L166_evidence.tsv`): derived at all nine L166-defining SNPs with zero
ancestral reads, and ancestral at all `G-Z6494`-defining SNPs with zero derived reads —
including the transversion `Z6494/FGC5674` (chrY:17131187 G>T) at 10 ancestral / 0 derived.

Furtwängler et al. 2020, *Ancient genomes reveal social and genetic structure of Late
Neolithic Switzerland*, Nat Commun 11:1915 (PMC7171184; author correction
10.1038/s41467-020-18561-y) reports Y assignments in Supplementary Table 1, made by visual
inspection of 1240k-captured SNPs against ISOGG v14.218. Relevant rows:

| site | n | terminal derived mutation | ISOGG longhand |
|---|---|---|---|
| Oberbipp Horgen | **7** | `PF3239` | G2a2a1a2a1 |
| Oberbipp Horgen | 2 | `FGC7739/Z6488` | G2a2a1a2a |
| Oberbipp Horgen | 1 | `L91/PF3246/S285` | G2a2a1a2 |
| Oberbipp Horgen | 2 | `PF3147` | G2a2a |
| Rapperswil | **1** (`SX10`) | **`L166`** | G2a2a1a2a1a |

`PF3239` is the direct parent of `L166` on the Iceman's own backbone
(`P15 > M3308 > PF3147 > PF3148 > PF3177 > L91 > Z6488 > PF3239 > L166`). A terminal-`PF3239`
call therefore means "derived at PF3239, nothing below demonstrated" — it is not an L166
exclusion. `SX10` establishes that `L166` is reachable in this dataset's panel and coverage
regime at least once.

A claim in public circulation states that 7 of 10 Oberbipp individuals are `G-L166`. The
published table does not say this; those seven are terminal `PF3239`. Third-party aggregators
relabel these samples downward (e.g. `MX211`, published as `PF3239`, is presented elsewhere as
`G-L166`). Whether the relabelling is nonetheless correct is an empirical question that the
source publication does not answer, and is the question this analysis addresses.

## Hypotheses

**H1 (derived).** The seven terminal-`PF3239` Oberbipp males are derived at the L166-defining
SNP set. The published stop at `PF3239` would then be coverage-limited rather than
phylogenetic, Oberbipp would be an L166-bearing population contemporary with the Iceman, and
the circulating claim would be accidentally correct.

**H2 (ancestral).** They are ancestral at the L166-defining SNP set. They sit on a sibling
branch under `PF3239`, the published assignment is exact, and the circulating claim is wrong.

**H0 (no power).** The L166-defining sites are uncovered or fail quality filters. The two
hypotheses are not distinguished.

H0 is a permitted and expected outcome. It is to be reported as H0, **not** as H2. Absence of
derived calls is not an ancestral call. This distinction is the single most important
commitment in this document.

## Predictions registered in advance

1. **Primary test is expected to be feasible.** `SX10` carries a published `L166` call from
   this same panel, so at least some L166-defining sites are 1240k targets. The seven Oberbipp
   males carry 50–391 Mbp each, more than `SX10`'s 79 Mbp in six of seven cases.

2. **The novel-variant test is predicted to fail for lack of power, in advance.** The Iceman's
   eight usable novel candidates (`results/iceman_y_novel_branch_candidates.tsv`) are
   uncatalogued by construction, therefore not 1240k targets, therefore reachable only by
   off-target reads. `SX10` has 79 Mbp of total sequence against `CGG017683`'s 1.18 Gbp, which
   itself yielded only 3 informative sites of 8 at 0.192x chrY. Expected informative sites here:
   **0–1**. If the result is zero, that is the predicted outcome and constitutes no evidence
   either way. It will not be reported as a negative finding.

3. **No prediction is registered on the direction of H1 vs H2.** The uniformity of the seven
   `PF3239` stops is equally consistent with a panel/coverage floor and with a real branch point.

## Data

Study `PRJNA608699` (ENA mirror), Furtwängler et al. 2020. All runs are
`Targeted-Capture` / `Hybrid Selection` / `SINGLE`, ~56 bp, UDG-half treated.
Raw FASTQ, not depositor-filtered BAMs.

| sample | run | reads | bases | role |
|---|---|---|---|---|
| MX210 | SRR11179329 | 7,105,270 | 391 Mbp | Oberbipp, terminal PF3239 |
| MX212 | SRR11179327 | 4,400,114 | 219 Mbp | Oberbipp, terminal PF3239 |
| MX211 | SRR11179328 | 3,455,110 | 186 Mbp | Oberbipp, terminal PF3239 |
| MX187 | SRR11179346 | 2,001,954 | 110 Mbp | Oberbipp, terminal PF3239 |
| MX182 | SRR11179349 | 1,673,721 | 85 Mbp | Oberbipp, terminal PF3239 |
| MX213 | SRR11179326 | 1,185,915 | 63 Mbp | Oberbipp, terminal PF3239 |
| MX209 | SRR11179330 | 876,635 | 50 Mbp | Oberbipp, terminal PF3239 |
| SX10 | SRR11179288 | 1,402,963 | 79 Mbp | Rapperswil, published terminal L166 |
| SX8 | SRR11179274 | 547,001 | 29 Mbp | mother of SX10; female, negative control |

Staged at `/mnt/AncientDNA/SwissLN-2020/`, MD5-verified against the ENA manifest.

## Marker sets

All coordinates GRCh38, from `results/iceman_y_L166_evidence.tsv`.

**Positive control — `PF3239` block.** The publication states these samples are derived here.
Recovering that independently validates download, alignment and pileup against a published
answer before any novel claim is made. A sample that fails this control is excluded from
interpretation.

**Primary — L166-defining set (9).** Transversions carry the evidential weight; the four
transitions are deamination-prone and are reported but not relied on alone.

| marker | chrY hg38 | anc>der | class |
|---|---|---|---|
| FGC5696 | 8525805 | C>A | transversion |
| Z6208 | 13776249 | G>A | transition |
| Z6219 | 13782251 | C>T | transition |
| Z6287 | 14368999 | C>T | transition |
| S19530/Z6213 | 15208473 | T>C | transition |
| FGC5721 | 16280455 | T>G | transversion |
| Z6516/FGC5675 | 19013205 | T>A | transversion |
| L166 | 21843737 | C>A | transversion |
| L167 | 21843756 | T>A | transversion |

**Exclusion set — `G-Z6494`-defining (3).** `FGC5687` 6019771 A>C, `Z6215` 11965520 T>C,
`Z6494/FGC5674` 17131187 G>T. The Iceman is ancestral at all three. Any Swiss sample derived
here is off his branch in the other direction and is reported as such.

**Secondary — Iceman novel candidates (8 usable).** chrY:7885869 T>A (transversion),
7899558 A>G, 8808031 C>T, 10768171 G>A, 10964462 T>C, 11414525 T>C, 11667647 T>C,
19647870 C>T. Marginal and MAPQ-rejected sites from the same table are excluded.

## Methods

Reference held constant with the Iceman analysis: `mapping/index/hg38p14DH3630O.fa` —
GRCh38.p14 + hs38d1 decoy + IPD-IMGT/HLA 3.63.0 (ID naming), 3.43 Gbp over 29,247 contigs,
chrY 57,227,415 bp. Verified to contain **zero oral/microbial contigs**.

The oral eHOMD decoy is deliberately **not** used. These are skeletal samples with no saliva
component, so the decoy has no contamination to absorb, and `METHODS.md` documents the cost of
including it anyway: Y-chromosome fragments absent from GRCh38 mismap readily to bacterial
genomes, so decoys can remove real chrY signal rather than only noise. Note the reference
naming is easy to misread — per `GRCh38_bwa_index.sh` the oral tag is emitted as
`O$(VERSION_ORAL | tr -d '.')`, so `...O915` carries eHOMD 9.15 while a bare trailing `...O`
carries none. Confirm composition from the `.fai`, not from the filename.

Alignment **differs by design** from the Iceman path and this is a declared deviation, not an
oversight. The Iceman is paired-end shotgun processed through the de-EAGER route
(`revert-bam.sh`, `util/eager_repair_bam.py`); these are 56 bp single-end capture reads with
no pairs to reconstruct. `bwa mem` is designed for reads ≥70 bp and loses sensitivity below
that, so `bwa aln`/`samse` is used instead, at the nf-core/eager 2.5.0 defaults
`-n 0.01 -k 2 -l 1024`. `-l 1024` exceeds any read length in this dataset and therefore
disables seeding, so the full read including its damage-bearing termini participates in the
alignment. Actual invocations are recorded in `RUNLOG.md` at run time.

Comparability between datasets is maintained downstream of alignment, not at it: identical
reference, identical marker sets, identical pileup thresholds. Uniformity of the alignment step
across incompatible library types would be cosmetic rather than methodological.

**Deduplication is mandatory**, via `samtools markdup -r`. This is capture data and therefore
PCR-amplified. The decision rules below count *independent* reads; without deduplication,
several duplicate copies of a single molecule would present as independent support and the
"≥2 reads" threshold would mean nothing. Duplication rates are recorded per sample.

Genotyping uses `annotate/y_markers_pileup.py --bam --ref --marker-file --label` for catalogued
markers, at the same settings as the Iceman run: `-q 25 -Q 20 -d 1000 --no-BAQ`.

The novel candidates have no catalogue names — being uncatalogued is what makes them novel — so
they are genotyped by position with `annotate/y_sites_pileup.py` at identical thresholds. That
tool also emits the mandatory `pct_mq0` / `n_mq60` audit columns and encodes the decision rules
below directly, so a call cannot be made by a different standard than the one registered here.

Mapping is performed by `mapping/map_se_adna.sh`, driven across the sample set by
`mapping/map_se_batch.sh`, and staging by `annotate/fetch_ena_runs.sh`. Genotyping across all
samples and marker sets is driven by `annotate/y_genotype_batch.sh`. The marker and site inputs
are committed under `markers/`. Every input to this analysis is therefore in the repository
rather than in transcribed shell history.

**Tool validation, performed before any Swiss read was mapped.** `annotate/y_sites_pileup.py`
was run against the Iceman BAM at all 21 sites in `markers/iceman_novel_candidates_all21.tsv`
and reproduced `results/iceman_y_novel_branch_candidates.tsv` — **21/21 verdicts**, including
the exact REJECT/MARGINAL classification. Output kept as
`results/iceman_y_novel_candidates_regen.tsv`. Two discrepancies against the older hand-built
table were found, both errors in that table rather than in the tool, and are recorded in
`NOTES.md`. Had the tool not reproduced those numbers, that would itself have been the finding
and this analysis would not have proceeded.

`annotate/y_genotype_batch.sh` was validated the same way and on the same schedule, against the
Iceman BAM: 9/9 `L166_defining` DERIVED with zero ancestral reads, 3/3 `Z6494_exclusion`
ancestral with zero derived reads, 10/10 `backbone_control` DERIVED. The entire published
G-L166-terminal finding is therefore recovered end to end by the committed tooling before this
dataset is touched, so a Swiss call is read against a chain known to work.

## Decision rules

Fixed in advance:

- **Derived call:** ≥2 independent reads supporting the derived allele with no conflicting
  reads; **or** 1 read at a transversion site with no conflicting reads, reported explicitly
  with its depth. A single transition read is never sufficient on its own.
- **Ancestral call:** ≥2 independent reads supporting the ancestral allele. A site with 1
  ancestral read is reported as `low_power`, not as ancestral.
- **Single-read C>T / G>A calls are treated as deamination artifacts, not evidence.** This
  failure mode is already documented in this repository for `CGG017682`, whose entire apparent
  G-panel signal was single-read transitions.
- **MAPQ audit is mandatory.** `pct_mq0` and `n_mq60` are reported per site. Depth alone does
  not qualify a site; 10 of 21 Iceman novel candidates were rejected on MAPQ at unremarkable
  depths, including two of three transversions previously presented as strongest. The threshold
  is **`pct_mq0` ≥ 30% → `REJECT_mapq`**, with `n_mq60 == 0` and single-strand support as
  `MARGINAL`. This cut was fixed during the Iceman analysis — before any read of this dataset
  existed in the working directory — where usable sites topped out at 14% MQ0 and rejected ones
  began at 39%; any cut inside that gap yields the identical split, and 30% is the round number
  in it. It is now encoded in `annotate/y_sites_pileup.py` rather than applied by hand.
- **Known no-go regions** are flagged, not silently included: chrY:11.1–11.7 Mb,
  chrY:56.69–56.88 Mb, ~26.6 Mb. These are emitted in a separate `region` column and are
  advisory: they do not by themselves reject a site, but no claim rests on a flagged site alone.
  Note that two of the eight "usable" Iceman candidates (`11414525`, `11667647`) fall inside the
  11.1–11.7 Mb window and were not flagged as such by the earlier hand-built table.
- **Output is per-SNP derived/ancestral/nocall.** No terminal haplogroup string is assigned to
  any sample. Assigning one is the error class this project exists to document.
- **Coverage is always reported with its denominator.** chrY depth on the ~12.57 Mb callable/MSY
  denominator, stated as such; the 59.37 Mb full-contig figure understates real site coverage by
  ~4.7x.

## Declared limitations

- **No pooled or cross-sample test is registered, and none will be run to rescue power.**
  Combining reads across the seven Oberbipp males to ask "are they L166" assumes they share a
  branch, which is the hypothesis under test. Pooling can support exactly one weaker claim —
  *is L166-derived present anywhere in the assemblage* — which a single clean derived
  transversion read in any one individual would establish. It cannot yield a per-individual
  assignment, a proportion, or a "dominated by" statement, and the circulating claim ("7 of 10")
  is per-individual, so only per-individual calls address it. If the per-sample tests come back
  H0, the answer is H0; pooling is not a fallback.
- **The seven may not be seven independent observations.** They are a collective cist burial with
  homogeneous Y assignments, which is the expected signature of a patrilineal kin group. If they
  are patrilineally related their Y chromosomes are shared by descent, the effective sample size
  is smaller than seven, and "7 of 10 carry the clade" reduces to "one lineage sampled seven
  times".

  **Resolved 2026-07-25 from Supplementary Table 5 — they are related.** `MX187`–`MX212` are
  first-degree and both are in this test set; `MX187`–`MX209` are second-degree. Further
  first-degree links tie the set to Oberbipp males not sampled here (`MX150`–`MX187`,
  `MX183`–`MX211`, `MX209`–`MX219`). Male–male first-degree pairs share a Y chromosome by
  descent. **No result here may be quoted as an independent observation, and no proportion may be
  computed over the seven.** Verdicts stay per-sample and are reported as what they are: calls on
  a kin group.

  A second consequence of that table matters more for interpretation. It carries a `same Y HG`
  column whose values include `same clade`, used where related males carry *different* published
  terminal labels — `MX150` (`L91`) with `MX187` (`PF3239`), `MX209` (`PF3239`) with `MX219`
  (`PF3147`). Men who share a Y by descent cannot differ phylogenetically, so those one-to-three
  node differences are depth, not branching. This is published evidence that the `PF3239` stop is
  coverage-limited — support for the *premise* of H1, obtained before any read of this dataset was
  called. It is **not** support for H1 itself: a floor being a floor says nothing about where the
  true terminal lies, which is exactly what the marker sets test. Registered here so it cannot
  later be presented as though the read-level data had supplied it.
- **"Same clade" is not "relative", and the relatedness claim is out of reach here.** Sharing
  `L166` implies a common patrilineal ancestor at or above the `L166` node, potentially a
  millennium or more before the Iceman. He is `L166*` — ancestral at all 93 markers in the
  subtree — so an L166-derived Oberbipp male could sit on any of the eight known sub-branches,
  i.e. a branch he is demonstrably not on. The test that would actually evidence recent shared
  ancestry is the shared-novel-variant test, already registered above as expected to yield 0–1
  informative sites. This analysis can therefore speak to clade membership and cannot speak to
  relatedness, and no result from it should be phrased as the latter.
- **The decision rules are asymmetric in sensitivity, and at this dataset's depth that biases
  which hypothesis is reachable.** A single transversion read calls DERIVED; an ancestral call
  needs ≥2 reads. So at DP 1 — which the first mapped samples show is the modal informative
  depth here, `MX182` at 0.0389x and `MX187` at 0.0545x chrY — H1 is detectable where H2 is not.
  If the truth is H1, thin coverage will still show it; if the truth is H2, thin coverage yields
  H0. **H0 must therefore not be read as leaning toward H2**, even though that is the direction
  the missing calls will appear to point. The asymmetry is deliberate and is not being changed:
  it exists because a single ancestral read is indistinguishable from a site that simply failed,
  whereas a transversion is damage-immune. Recorded here, before the results, so that it cannot
  later be mistaken for a post-hoc explanation of a null. Registered 2026-07-25 with two of nine
  samples mapped and no L166-defining site yet called in any of them.
- 56 bp single-end capture; no insert-size model, no mate rescue.
- `bwa aln` is not ALT-aware, so `bwa-postalt.js` is not applied and the 23 MB `.alt` companion
  to the reference goes unused. Accepted for chrY marker genotyping, where MAPQ filtering does
  the equivalent work; recorded because it is a real difference from the Iceman processing.
- 1240k capture means off-target coverage is the only route to uncatalogued sites.
- UDG-half treatment leaves residual damage at read termini.
- Sample sizes are single individuals. Nothing here supports a population-level claim.

## Out of scope

This analysis tests Y-chromosome marker status. It does not, and cannot, establish where the
Iceman came from or what culture he belonged to. Shared terminal Y status at this depth is
compatible with a common paternal ancestor centuries or millennia earlier. Origin questions are
addressed by isotopic and autosomal evidence — for the Iceman, Müller et al. 2003, *Science*
302:862, restricts his childhood range to valleys within ~60 km south(east) of the discovery
site. No result from this analysis will be presented as bearing on that.
