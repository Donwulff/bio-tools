# Pre-registration: read-level status of the unhedged `G-L166` individuals

Written 2026-07-26, **before any read from any of these five samples has been staged, mapped or
examined**. Hypotheses, decision rules, power predictions and declared limitations below are
fixed as of that date. Results go in `RUNLOG.md` and `NOTES.md`; this file is not to be edited
to match them.

This is a **separate document rather than a fourth amendment** to
`PREREG_swiss_neolithic_L166.md` because the dataset has left Switzerland: the samples are
Hungarian, Sardinian and Bavarian, from three unrelated publications. It **inherits by
reference** and does not restate: the decision rules, the MAPQ audit and its 30% cut, the no-go
regions, the per-SNP output rule, the coverage-denominator rule, and the standing prohibition on
cross-sample pooling. Where this document is silent, that one governs.

## Background and what is actually being tested

The Iceman is terminal `G-L166` by read-level evidence. Work to date has tested the *hedged*
label: the compilation at haplogroup.info writes ISOGG `G2a2a1a2a1(a)` — parenthesised, meaning
provisional — for 21 Swiss Late Neolithic males, and a claim in circulation reads seven of the
Oberbipp subset as `G-L166`. Fifteen Oberbipp-area samples have now been genotyped. Across four
mutually independent Y lineages, the two defining transversions returned **21 reads, all
ancestral, none derived** (`results/swiss15/`, `NOTES.md` 2026-07-26). The hedged label does not
verify.

That leaves the complementary question, and it is the one with the sharper edge. A smaller set
carries the **unhedged** longhand `G2a2a1a2a1a` with bare `G-L166` in both the YFull and
Y-Haplotree columns — no parentheses, no asterisk. If the hedging convention carries the meaning
this project has inferred for it, these should verify at read level where the hedged ones did
not. **The convention itself is therefore under test here, not just the samples.** That
inference has never been confirmed against a published legend; no legend has been found. This is
the first test that could support or break it.

### The set is six individuals, not seven

Enumerating the compilation for unhedged `G2a2a1a2a1a` / `G-L166` returns seven rows. Two of
them are the same man:

| field | `UNTA58_68Sk1` | `E09538` |
|---|---|---|
| ¹⁴C lab number | **MAMS-29075**, 3870±30 BP | **MAMS-29075**, 3870±30 BP |
| calibrated | 2462–2209 calBCE | 2464–2210 calBCE |
| mtDNA | J1c | J1c |
| coordinates | 48.31611, 10.892 | 48.32, 10.89 |
| colloquial | `UNTA58_68Sk1` | "Feature 68 Skeleton 1" |
| site | Haunstetten – Unterer Talweg 58–62 | Unterer Talweg 58-62, Augsburg, Bavaria |
| source | MittnikScience2019 | OlaldeNature2018 |

A radiocarbon lab number is assigned to one bone sample, and the Mittnik identifier decodes to
the Olalde description (site UNTA58, feature 68, skeleton 1). Haunstetten is a district of
Augsburg. These are two publications of one individual under two aggregation IDs.

Registered here because it is a **counting** finding, fixed before any read: any tally over this
list — including any this project publishes — is otherwise one too high. `SX10` (`TU876`) is the
seventh row and has already been genotyped, leaving **five new individuals**.

## Hypotheses

**H1 (verified).** The unhedged samples are derived at the L166-defining set with no conflicting
ancestral reads. The published unhedged label is exact, the hedged/unhedged distinction is
meaningful, and genuine `G-L166` individuals exist outside the Iceman — in Hungary at
near-contemporary date, in Chalcolithic Sardinia, and in Bavarian Bell Beaker.

**H2 (falsified).** They are ancestral at the L166-defining set. The unhedged label is
over-called in the same way the hedged one is, the parenthetical convention carries no
read-level meaning, and the whole `G-L166` attribution layer in this compilation is unreliable
rather than merely provisional in places.

**H3 (split).** Some verify and some do not. This is a distinct outcome and must not be
collapsed into either of the above; it would mean the label's reliability is sample-specific and
would require per-sample reporting with no summary verdict.

**H0 (no power).** The defining sites are uncovered or fail quality filters, and the hypotheses
are not distinguished.

H0 remains a permitted, first-class outcome and is to be reported as H0, never as H2. Unlike the
Swiss test, however, H0 is **not** the expected outcome here — see the power prediction below,
which is registered precisely so that an H0 result counts as a surprise requiring explanation
rather than as the anticipated null.

## Power prediction, registered in advance

This is the material difference from the Swiss analysis and it is settled **before staging**, not
discovered afterwards.

Four of the five samples are targeted-capture libraries. Autosomal coverage is irrelevant to
whether a specific site was interrogated: on capture data, depth at an off-panel position is
whatever off-target carryover delivers, so a 3.9x sample can still carry zero reads at the one
transversion the test turns on. Panel membership was therefore resolved first, from the panel
definition itself, by `annotate/panel_membership.py` — committed tooling, not a hand lookup —
with the result kept as `results/panel/1240k_marker_membership.tsv`.

Panel: `51.2.2M.snp` (EIGENSTRAT, 2,144,502 sites, 32,681 on chrY; md5 `6389a9c8…`), lifted
GRCh38→GRCh37 from `hg38ToHg19.over.chain`.

**18 of our 22 markers are on-panel.** Decisive for this test:

| marker | class | on panel | panel identifier |
|---|---|---|---|
| `L166` | C>A **transversion** | **yes** | `L166` (23989884) |
| `L167` | T>A **transversion** | **yes** | `L167` (23989903) |
| `Z6219` | C>T transition | **yes** | `snp_24_15894131` |
| `PF3239` | C>T transition | **yes** | `snp_24_17317628` |
| `Z6494` | G>T transversion | **no** | — |
| `Z6208`, `Z6488`, `FGC5687` | — | **no** | — |

Matching was done by lifted position, not by name. Name matching alone returns 5 of 22 and is
**wrong**: panels carry synonyms, and `P15`, `M3308`, `L91` and `P287` appear as `PF3112`,
`L1259`, `S285` and `rs4116820` respectively. All 18 positional hits agree on both alleles,
which is the check that the liftover is correct rather than coincidental.

Three consequences are registered as predictions:

1. **The primary test is expected to have real power.** Both defining transversions are
   deliberate capture targets, and total sequence per sample (167–849 Mbp) equals or exceeds the
   best Oberbipp library (`MX210`, 391 Mbp). **I expect calls, not H0.** If the primary sites
   come back uncovered in samples of this depth with the sites on-panel, that is itself a finding
   about the libraries and will be reported as such.

2. **The `Z6494` exclusion test is degraded and this is stated now, not later.** Its only
   transversion, `Z6494` itself, is off-panel, as is `FGC5687`. On these libraries the exclusion
   set reduces to `Z6215` (T>C, a **transition**, deamination-exposed). A `G-Z6494` exclusion
   therefore **cannot be made to the standard applied to the Iceman** in these samples. Where the
   exclusion is untestable it will be reported as untestable, and no sample will be called
   `L166*` — terminal, no known subclade — on this data. The Iceman's `L166*` status rests on his
   own 10-read ancestral `Z6494` transversion and is untouched by this limitation.

3. **`Z6208` off-panel means the `G-Z6208` question is out of reach here.** Since the published
   Iceman label `G2a2a1a2a1a1b (G-Z6208*)` is the very claim this project corrects, note that
   these samples cannot speak to it in either direction.

## The `Z6219` conflict — what this dataset can and cannot do

`PREREG_swiss_neolithic_L166.md` Amendment 2 §C left `Z6219` unresolved: derived in every
Oberbipp male while `L166` was ancestral in the same individuals, with two live explanations —
**(i)** the catalogue mis-levels `Z6219`, which actually sits at or above `PF3239`; **(ii)**
`Z6219` marks a real node between `Z6488`/`PF3239` and `L166`. Test C returned H0 exactly as its
registered prior predicted, because `MX204` — the only admissible candidate — had no reads at
`PF3239` and so failed the admissibility gate.

**These five samples cannot settle it either, and that is registered now so no post-hoc reading
of their data can claim otherwise.** Explanations (i) and (ii) make the *same* prediction for a
genuine `L166` man: derived at both `Z6219` and `L166`. Discriminating them requires an outgroup
**ancestral at `PF3239`**, and every sample here is published at or below `L166`. No such
outgroup is present.

Two things are nonetheless registered as decision-relevant:

- **A refuting outcome exists and is defined in advance.** A sample called **derived at `L166`
  and ancestral at `Z6219`** (≥2 ancestral reads, 0 derived) would refute *both* (i) and (ii),
  placing `Z6219` on a sibling branch rather than anywhere above `L166`. It is not expected. It
  is registered so that it would be recognised rather than explained away.
- **`PF3239` being on-panel opens the route that Oberbipp lacked.** The admissibility gate that
  `MX204` failed is a deliberate capture target on this panel. This does not help within the
  present five, but it means the outgroup test is tractable on *any* 1240k-captured G2a male
  ancestral at `PF3239`. That search is **explicitly out of scope for this pre-registration** and
  requires its own, written before those samples are selected — selecting outgroups after seeing
  which ones would settle the question is the failure mode this discipline exists to prevent.

## Data

Runs identified from ENA by `sample_alias`; **no reads staged at the time of writing.** Study
accessions were verified by query, not assumed — an earlier accession in `RUNLOG.md`
(`PRJEB36959`) proved to be a *Leptospira* 16S study, and two guesses made while locating these
five resolved to a *Medicago* rhizosphere GWAS and a feline cardiomyopathy study.

| sample | study | source | runs | total bases | autosomal | site, date |
|---|---|---|---|---|---|---|
| `I5118` | `PRJEB23635` | Olalde 2018 | `ERR2207344`, `ERR2207549` | 375.5 Mbp | 1.631 | Mezőcsát-Hörcsögös, HU; 3300–3000 BCE |
| `I14678` | `PRJEB35980` | Fernandes 2020 | `ERR3800865`, `ERR3800866` | 848.8 Mbp | 3.867 | Serra Crabiles, Sardinia; 2454–2201 calBCE |
| `I14677` | `PRJEB35980` | Fernandes 2020 | `ERR3800863`, `ERR3800864` | 717.2 Mbp | 3.204 | Serra Crabiles, Sardinia; 2464–2294 calBCE |
| `I15942` | `PRJEB35980` | Fernandes 2020 | `ERR3800873`, `ERR3800874` | 315.0 Mbp | 0.957 | Anghelu Ruju, Sardinia; 2459–2209 calBCE |
| `UNTA58_68Sk1` | `PRJEB34400` | Mittnik 2019 | `ERR3518170`, `ERR3518171` | 166.6 Mbp | 1.040 | Augsburg, Bavaria; 2462–2209 calBCE |

All single-end. `I5118` is deposited as `library_strategy=OTHER`; the three Fernandes runs as
`Targeted-Capture`/`other`; the Mittnik runs as `Targeted-Capture` with
`library_selection=RANDOM`, which is internally contradictory metadata and will be resolved from
the observed on/off-panel depth ratio rather than from the deposited annotation.

**`I5118` is the sample of primary interest.** At 3300–3000 BCE it is near-contemporary with the
Iceman (~3350–3100 BCE) and is the only non-Sardinian, non-Bell-Beaker individual in the set.

## Independence

`I14677` and `I14678` are from the **same tomb** — Serra Crabiles, Tomb III Cella A, identical
coordinates — and both carry `Kinship-Notes = n/a`. That field records *absence of a documented
relationship*, which is not evidence of unrelatedness, and collective tomb burial is precisely
the context in which patrilineal kin cluster. Their mtDNA differs (`K1a1b1` vs `V`), which
excludes a maternal sibship but says nothing about the Y.

**Registered rule:** `I14677` and `I14678` are reported per-sample but are **not** counted as two
independent Y observations unless their reads demonstrate it. If both verify at `L166`, the pair
contributes **one** observation to any count of independent lineages, and the possibility that
they are one patriline is stated wherever the count appears. This is the Oberbipp Family A lesson
applied before the fact rather than after it.

The remaining three are from different countries and unrelated publications and are independent.

## Marker sets and methods

Marker sets are unchanged and are the committed files `markers/L166_defining.txt`,
`markers/Z6494_exclusion.txt`, `markers/backbone_control.txt` — with the `Z6494` degradation
declared above. Reference, aligner, deduplication and pileup settings are inherited unchanged
from `PREREG_swiss_neolithic_L166.md` so that results are comparable to `results/swiss15/`:
`mapping/index/hg38p14DH3630O.fa`, `bwa aln`/`samse` at nf-core/eager 2.5.0 defaults,
mandatory `samtools markdup -r`, pileup at `-q 25 -Q 20 -d 1000 --no-BAQ`.

Staging by `annotate/fetch_ena_runs.sh`, mapping by `mapping/map_se_batch.sh`, genotyping by
`annotate/y_genotype_batch.sh`, verdicts by `annotate/y_prereg_verdict.py`. Outputs to
`results/unhedged/`. `results/swiss/` and `results/swiss15/` are not to be modified.

The backbone control (`P15 → M3308 → PF3147 → PF3148 → PF3177 → L91 → PF3239 → P287`, all
on-panel except `Z6488`) governs as before: a sample failing it is excluded from interpretation
rather than reported as a negative.

## Declared limitations

- **Read depth per sample is unknown at registration.** Total bases are known; chrY depth is not,
  and capture efficiency varies by library. The power prediction above is a prediction, not a
  measurement.
- **No `L166*` (terminal) call is possible on these libraries**, for the `Z6494` reason given
  above. The most this analysis can establish for any sample is `L166`-derived.
- **Verifying the label does not establish relatedness to the Iceman.** Shared `L166` implies a
  common patrilineal ancestor at or above that node, potentially a millennium or more before any
  of these individuals. The Iceman is `L166*` and an `L166`-derived sample here could sit on any
  of the eight known sub-branches — branches he is demonstrably not on. This constraint is carried
  over verbatim in force from the Swiss document and applies identically.
- **Five individuals from three studies support no population-level claim** about the
  distribution of `G-L166` in Chalcolithic Europe, and none will be made. Absence of `G-L166` in
  any region is not tested by this design and cannot be inferred from it.
- **The hedged/unhedged convention remains an inference.** If these five verify, that is
  *consistent* with the convention meaning what this project reads it as meaning; it is not proof,
  because no published legend has been located. The convention's status will continue to be
  reported as inferred.
- The asymmetry of the decision rules (1 transversion read → DERIVED; ancestral requires ≥2) is
  unchanged and still means H1 is more reachable than H2 at low depth. At the depths expected here
  it should bite less than it did at Oberbipp, but **H0 must still not be read as leaning toward
  H2**.

## Out of scope

Unchanged from `PREREG_swiss_neolithic_L166.md`: this tests Y-chromosome marker status and cannot
establish origin, migration or cultural affiliation for any individual, including the Iceman. The
geographic spread of this set — Hungary, Sardinia, Bavaria — makes that temptation larger, not
smaller, and no result here will be presented as bearing on where anyone came from.
