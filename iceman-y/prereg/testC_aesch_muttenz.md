# Pre-registration: Test C — does the `Z6219` node recur outside Oberbipp?

Written 2026-07-26, **before any Aesch or Muttenz read has been staged, mapped or genotyped**.
This is the "own power statement" that `iceman-y/prereg/Z6219_node.md` §7 requires before Test C may run:

> **Test C — new samples.** Out of scope for this registration. The 21 hedged Aesch/Muttenz
> individuals are named in P3 as the natural cohort, but selecting and staging them is a separate
> exercise requiring its own power statement, and nothing in this document authorises it.

It is a new registration rather than an extension because it introduces a cohort with its own power
question, which is the third row of the table in `iceman-y/PROTOCOL_extending_analyses.md`.

---

## 1. Why this test exists

The registered verdict of `iceman-y/prereg/Z6219_node.md` (2026-07-26) is **H1 supported**, recorded together
with the limit that makes it fragile:

> **One site.** All four "independent lineages" are Oberbipp. H3 predicts the pattern is confined to
> that lineage; H1 predicts it recurs elsewhere. **No second population has been tested.** E2
> attempted one (`CGG017683`) and returned H0 on depth. This is the largest single weakness.

H3 — homoplasy, or a variant private to one patriline — survives everything done so far, because
everything done so far is one patriline at one site. Test C is the only registered test that
separates H1 from H3. That is its entire purpose, and if it cannot do that it has failed regardless
of what else it produces.

## 2. The cohort, and three corrections it forced

**The phrase "21 hedged Aesch/Muttenz individuals", used in `iceman-y/prereg/Z6219_node.md` §7 and in
`RUNLOG.md`, is wrong, and this document does not inherit it.** Enumerating the All Ancient DNA
compilation (`/mnt/AncientDNA/all-ancient-dna.2026-07-26.txt`, 16,972 rows) for
`Y-Haplotree-Variant == "G-L166*"` gives **21 individuals** — but they are:

| site | n | status |
|---|---|---|
| Aesch (CH) | 13 | untested |
| Muttenz (CH) | 2 | untested |
| Oberbipp Horgen (CH) | 6 | **already genotyped**, in `iceman-y/results/swiss15/` |

So the 21 is a real count, but 6 of them are the Oberbipp men this project has already tested, and
the new cohort is **15 individuals, not 21**. The error was reading "21 hedged" and "Aesch/Muttenz"
as describing the same set.

**Second correction.** `NOTES.md` (line 897) states that Aesch and Muttenz "do not appear at all" in
the compilation. That is true only of the `G-PF3239` query it was written under. Both sites are
present in quantity under the IDs `Aes1`–`Aes25` and `RA42`–`RA64`, which no search for the strings
"Aesch" or "Muttenz" in an ID column will match. This is the same class of error as the `G-L166`
literal-string search recorded in the 2026-07-26 correction immediately above it in that file.

**Third correction, and it constrains what this test may claim.** P3 describes the cohort as
"Published-`PF3239` males from Aesch and Muttenz". The cohort registered here is **not** that set.
It is defined by the *compilation's* `G-L166*` label, and the compilation's labels are a hobbyist
re-derivation, already documented in `NOTES.md` as disagreeing with the publication in both
directions. The compilation's own `Y-DNA` column is likewise re-derived, not published: it calls
`MX187` `G2a2a1a2a1a` where the paper says `PF3239`, and `MX183` `G2a2a1a2a1` where the paper says
`Z6488`. **No claim in this document rests on any sample being published `PF3239`**, and the
publication's Supplementary Table 1 has not been re-read for Aesch at the time of writing. See §7.

The cohort, its per-sample coverage proxy and its kin-group assignment are committed in
`iceman-y/results/testC_power/candidates.tsv` and `iceman-y/results/testC_power/candidates_lineage.tsv`, both
regenerable from the compilation file by the command in `RUNLOG.md`.

### The cohort is not 15 independent observations

Oberbipp's lesson was that seven men were largely one family. Aesch is the same shape, and the
kinship is documented in the compilation's own `Kinship-Notes`:

| kin group | candidates in it | documented structure |
|---|---|---|
| Aesch Family D (12 members) | `Aes1`, `Aes12`, `Aes19`, `Aes20`, `Aes21`, `Aes23` | `Aes12`–`Aes19` 1st degree; all others 2nd or 3rd degree |
| Aesch Family C (4) | `Aes14`, `Aes22` | 2nd or 3rd degree |
| Aesch Family A (3) | `Aes6`, `Aes17` | 2nd or 3rd degree |
| Aesch Family B (3) | `Aes4`, `Aes13` | 2nd or 3rd degree, via `Aes18` (female) |
| — | `Aes7` | no documented kinship |
| — | `RA61`, `RA62` | no documented kinship |

**Seven kin groups at two sites**, against Oberbipp's one site and — after Family A pooling — very
few. That, and not the raw individual count, is what Test C buys.

> **Annotation added after the §7 checks, 2026-07-26.** The table above is sourced to the
> compilation and it is *not* the publication's. Furtwängler 2020 Supplementary Table 5 documents
> only one of these links among the candidates (`Aes12`–`Aes19`). See "§7 checks — completed" at the
> end of this document. The registered pooling rule is unaffected; the attribution is.

## 3. Hypotheses

**H1 — `Z6219` marks a real node between `PF3239` and `L166`.** Then the `Z6219`-derived,
`L166`-ancestral state is a property of a lineage that was geographically widespread, and should be
found in kin groups at Aesch and/or Muttenz that are not the Oberbipp patriline.

**H3 — homoplasy, or an Oberbipp-private variant.** Then the `Z6219`-derived state does not recur.
Callable Aesch and Muttenz men are `Z6219`-ancestral, or `Z6219`-derived only where they also turn
out `L166`-derived (i.e. genuinely below the node, telling us nothing about an intermediate branch).

**H0 — no power.** Too few individuals reach the registered depth threshold to distinguish these.
**This is the single most likely outcome and is a first-class result**, not a failed run. E2 returned
H0 and is committed for that reason.

## 4. Power, fixed before any read is staged

Computed by `annotate/y_power_estimate.py`, committed before it was run against the candidates.
Outputs in `iceman-y/results/testC_power/`. The estimate is empirical throughout: no model of capture
efficiency, probe design or mappability appears anywhere in it. Whatever those contribute is already
inside the measured per-marker read rate, which is the point of measuring rather than modelling.

**Calibration.** Fifteen already-mapped Swiss libraries from the same study, the same capture
protocol and the same instrument (`Targeted-Capture` / `Hybrid Selection`, HiSeq 3000, verified in
the ENA metadata for both cohorts) yielded **135,398 chrY reads at MAPQ ≥ 25**, and against that
base:

| marker | reads across all 15 | rate per chrY read | 95% Poisson interval |
|---|---|---|---|
| `L166` | 14 | 1.03e-4 | 5.6e-5 – 1.7e-4 |
| `Z6219` | 10 | 7.39e-5 | 3.5e-5 – 1.4e-4 |
| `L167` | 7 | 5.17e-5 | 2.1e-5 – 1.1e-4 |
| `FGC5671` | 5 | 3.69e-5 | 1.2e-5 – 8.6e-5 |

**Only 5 of the 22 YFull `L166`-defining positions attracted a callable genotype in any of the 15
libraries**, and 10 of the 22 attracted zero reads in all 15 combined — including `FGC5696`,
`FGC5721`, `Z6516` and `S19530`, all of which *are* on the 1240k panel. Panel membership therefore
does not predict coverage at this locus, and the testable marker set is effectively `L166`, `L167`
and `Z6219`. Pooling more markers to buy power is not available.

**Coverage proxy.** A candidate's mapped chrY yield is unknown until it is mapped, which is the
expense this estimate exists to avoid, so it is predicted from the compilation's `NRY` column. Over
the 12 calibration samples carrying both, `chrY_reads_mq25 / NRY` runs **0.543 – 0.799** with median
**0.6514** — tight enough to predict on, and the ratio actually used.

One calibration row is corrupt: `MX182` carries `NRY = 2` against 9,148 mapped chrY reads, a defect
recorded in `NOTES.md` on 2026-07-25, before this test existed. A through-origin fit would let that
single row raise the ratio to 0.7455, a 14% optimistic bias, because a sample with a near-zero
denominator and real reads dominates the numerator. The tool therefore uses the **median of
per-sample ratios**, which cannot be moved that way. `MX182` is not dropped — it stays in the
committed table, visible, and simply stops dominating the fit. This is a robustness choice made on
the estimator, not a curation of the data, and it was made before any candidate power number was
produced.

**Self-check against the cohort it was fitted on.** Run against the Oberbipp/Rapperswil samples
themselves, the model predicts **2.08** callable at `Z6219` (observed: **2** — `MX210`, `SX10`) and
**3.15** at `L166` (observed: **4**). It reproduces the cohort it came from.

**Result for the 15 candidates**, at the registered depth threshold of ≥ 2 independent reads:

| quantity | expectation | P(≥1) | P(≥2) |
|---|---|---|---|
| individuals callable at `Z6219` | **4.0** | 0.994 | 0.953 |
| individuals callable at `L166` | **5.9** | 1.000 | 0.997 |
| individuals callable at **both** | **2.1** | 0.918 | 0.666 |

Callable-at-both is the number that matters, because the intermediate-branch signature is
`Z6219`-derived **and** `L166`-ancestral **in the same man**. One marker alone cannot place anyone
between two nodes.

**Power to observe the split P3 predicts**, as a function of the unknown derived frequency in this
population:

| derived frequency | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 |
|---|---|---|---|---|---|
| P(both alleles seen at `Z6219`) | 0.33 | 0.63 | **0.79** | 0.63 | 0.33 |

**The dominant uncertainty is the marker rate, not the cohort.** Re-running at the bounds of the
`Z6219` Poisson interval (`--rate-scale 0.479` and `1.838`) gives expected callable of **1.3** and
**7.9**, and split power at frequency 0.5 of **0.23** and **0.98**. The honest statement is that
Test C's power is somewhere between "barely worth running" and "decisive", and which one will not be
known until the reads are mapped. It is registered anyway, because the alternative is to keep the
verdict resting on one site.

## 5. Falsifiers

**FC1 — H1 falsified.** Two or more independent kin groups outside Oberbipp are callable at both
`Z6219` and `L166` and are **`Z6219`-ancestral** while `L166`-ancestral. A man ancestral at both is
below `PF3239` and above the `Z6219` node; enough of them, with no derived counterpart anywhere
outside Oberbipp, and the intermediate branch is an Oberbipp property, which is H3.

**FC2 — the node is not where it was placed.** Any individual **`Z6219`-ancestral and
`L166`-derived**, at depth, in either direction of the damage rules. That ordering is impossible
under H1 and would falsify the `PF3239 → Z6219 → L166` order outright rather than merely weakening
it. This is the cleanest falsifier in the document and it is symmetric: it is exactly as likely to
appear if I am wrong as if I am right.

**FC3 — the cohort is one lineage.** If every callable candidate turns out to carry an identical Y
haplotype across the backbone control set, the kin-group structure in §2 is not Y-independent and
the whole cohort collapses to one observation, exactly as Oberbipp did. Then the outcome is **H0**,
not support for anything.

## 6. Predictions if H1 is true

**PC1.** At least one kin group outside Oberbipp is `Z6219`-derived and `L166`-ancestral.

**PC2.** Callable candidates **split** at `Z6219`. Uniformity in either direction is uninformative
rather than confirmatory, as already registered in P3 — with ~4 callable individuals, all-ancestral
and all-derived are both entirely ordinary outcomes under H1.

**PC3.** The 10 YFull positions that attracted zero reads across 15 Oberbipp libraries attract
approximately zero here too. This is a boring prediction and it is registered because it is the
control: if they suddenly produce coverage, something is wrong with the mapping, not with the
haplotree.

## 7. Required checks before any genotype is read

Each of these can invalidate the cohort, and each is a provenance or staging check, which
`iceman-y/PROTOCOL_extending_analyses.md` explicitly exempts from registration. They are listed here so that
skipping one is visible.

1. **ID mapping.** The compilation uses `Aes12`; ENA uses `Aesch12`; the compilation uses `RA61`,
   ENA uses `SNPRA61`. Numbering is *assumed* to correspond and has not been verified. A silent
   off-by-one here would genotype the wrong individuals and produce a confident answer. Verify
   against Furtwängler 2020 Supplementary Table 1 before mapping.
2. **Sex.** `RA61` and `RA62` carry an empty `Sex` field in the compilation. Confirm both are male
   before treating a Y genotype from them as meaningful.
3. **Published labels.** Re-read Supplementary Table 1 for Aesch and Muttenz. §2 records that this
   cohort is *not* the published-`PF3239` set P3 named, and the actual overlap is unknown.
4. **`Aesch13` has `NRY = 17`** and predicted chrY coverage of ~11 reads. It is registered as a
   candidate for completeness and is expected to return no power. Report it, do not quietly drop it.
5. **Backbone control first.** `markers/backbone_control.txt` must recover before any
   `L166`-defining call from these BAMs is believed, exactly as for the Swiss 15.

## 8. Counting rules, fixed before the scan

- **Depth threshold is unchanged**: ≥ 2 independent reads for a call in either direction. `Z6219` is
  a C>T deamination-prone transition, so a single derived read is `nocall_damage_prone_1read` and
  **is not evidence**. `L166` and `L167` are transversions and damage-immune.
- **`MAPQ ≥ 25` stands**, and with it the consequence registered in `iceman-y/prereg/uniqueness_filter.md`:
  every `L166` derived count is a **lower bound**, and `no_coverage` at `L166` is not evidence of
  anything. Uniqueness filtering was tested and rejected on 2026-07-26; it is not revisited here.
- **Pooling.** Only documented **1st-degree male–male** pairs may be pooled, because only those
  share a Y chromosome by descent. Among candidates that is `Aes12`–`Aes19` and nothing else. A
  2nd- or 3rd-degree relative may be linked through a female and must not be pooled; pooling across
  an unshared Y would manufacture a false split inside one unit.
- **Independence.** Kin groups negative at 3rd degree are **not** thereby Y-independent — four
  families at one Neolithic village can share a patriline beyond the resolution of READ/lcMLkin.
  Counts of "independent lineages" in any write-up must carry that caveat. This is why FC3 exists.
- **Any decisive C>T or G>A call requires read-terminus and library damage-profile evidence**
  (`annotate/y_read_evidence.py`), per `iceman-y/prereg/Z6219_node.md` §8. That is protocol, not a rescue
  applied when a result depends on it. Note `--damage-profile` is a *mode switch* that suppresses
  the per-read dump, so it takes two invocations.
- **`untested` is not `negative`.** The failure mode named in `iceman-y/PROTOCOL_extending_analyses.md` §2.

## 9. Decision rules

Outcome to be recorded as exactly one of:

- **FC1 obtains → H1 weakened, H3 favoured.** The `H1 supported` verdict in `iceman-y/prereg/Z6219_node.md`
  is amended, not deleted, and the amendment names this test.
- **FC2 obtains → the node order is falsified.** This overrides everything else in the document.
- **PC1 met → H1 supported at a second site.** The single-site limitation recorded in the verdict is
  discharged to the extent of the number of independent kin groups actually observed, which must be
  stated as a number and not as "replicated".
- **H0 — no power.** Fewer than two candidates callable at both markers, or FC3 obtains. Recorded
  and committed, as E2 was.

No terminal haplogroup label is assigned to any individual by this analysis. Per-SNP output only.

## 10. Peeking declaration

**Read before writing this document**, and it is why the document exists: the compilation's
`Y-Haplotree-Variant`, `NRY`, `Location` and `Kinship-Notes` columns for all Aesch and Muttenz rows;
ENA read counts, library strategy and instrument for all 97 runs of `PRJNA608699`; and the per-marker
**depths** in `iceman-y/results/testB/swiss_yfull_L166_defining.tsv`. Depths were read to establish whether
the test is possible at all, which `iceman-y/PROTOCOL_extending_analyses.md` permits explicitly.

**The Oberbipp and Rapperswil allele calls at `Z6219` and `L166` have been read** — they are the
result this test exists to check, quoted throughout `NOTES.md`, and were read long before this
document.

**Not read, and not in existence at time of writing:** any Aesch or Muttenz sequence data. No run
has been fetched, no BAM built, no genotype called. The power numbers in §4 were produced entirely
from read counts, an external coverage proxy and the Oberbipp marker rates.

## 11. What this analysis will not claim

- Nothing about the Iceman's own placement. `iceman-y/prereg/Z6219_node.md` §10 governs and is unchanged.
- Nothing about the Horgen culture, or about any site's cultural attribution. The compilation's
  `Culture_Grouping` assigns "Horgen culture" to Aesch, Muttenz, Oberbipp and Rapperswil alike, which
  the publication does not support, and this document uses it for nothing.
- **Nothing about whether the compilation's `G-L166*` labels are correct.** They are the *sampling
  frame*, chosen because it is reproducible from a committed file. Confirming or refuting a hobbyist
  re-derivation is not a finding about the haplotree.
- No re-interpretation of any prior finding in the same commit as the results. Calls are recorded
  first; what they mean is a separate question asked afterwards.

---

## §7 checks — completed 2026-07-26, before any read was staged

Run against Furtwängler 2020's corrected Supplementary Information PDF, fetched and read directly
(`static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15560-x/MediaObjects/41467_2020_15560_MOESM1_ESM.pdf`,
25 pp., Supplementary Tables 1 and 5). Not from a search summary and not from the compilation.

**The cohort survives all three checks. Nothing in §§1–11 above is amended.**

### 1. ID mapping — confirmed, and it runs the opposite way from the assumption

Supplementary Table 1 uses `Aesch1`, `Aesch12`, `Aesch20`, `RA61`, `RA62` — **identical to the ENA
`sample_alias` values**. It is the *compilation* that abbreviates to `Aes1`, `Aes12`, `RA61`. So the
mapping applied when building `candidates.tsv` (`Aes` → `Aesch`, `RA` → `SNPRA`) recovers the
publication's own IDs, and the ENA↔publication correspondence needs no inference at all. The
off-by-one hazard named in §7.1 does not exist for this cohort.

### 2. Sex of `RA61` and `RA62` — confirmed male

Supplementary Table 1 is titled *"Y chromosomal haplogroup assignment for **all male individuals**"*.
Both appear in it. The empty `Sex` field in the compilation is a gap in the compilation, not an
ambiguity in the source.

### 3. Published labels — the overlap with P3's set is now known, and the cohort is a superset

Of the 15 registered candidates, **12 are published `PF3239`** and 3 are published *shallower*:

| candidate | published terminal | published YHG |
|---|---|---|
| `Aesch6` | `PF3147` | `G2a2a` |
| `Aesch7` | `FGC7739/Z6488` | `G2a2a1a2a` |
| `Aesch20` | `FGC7739/Z6488` | `G2a2a1a2a` |

The other 12 — `Aesch1`, `Aesch4`, `Aesch12`, `Aesch13`, `Aesch14`, `Aesch17`, `Aesch19`,
`Aesch21`, `Aesch22`, `Aesch23`, `RA61`, `RA62` — are published `PF3239` / `G2a2a1a2a1`.

That is **all 10 published-`PF3239` individuals at Aesch and both at Muttenz**, confirming the
per-site counts recorded in `NOTES.md`. The registered cohort is therefore a strict **superset** of
the set P3 named, and Test C can report under either definition without re-registration. The three
shallower individuals are retained: they were selected by a stated, reproducible criterion, and
dropping them now on a label read afterwards is rule-shopping. `Aesch6`'s predicted power is 0.036
at `Z6219`, so retaining it costs nothing either way.

### 4. Kinship — the registered pooling rule is confirmed by the source, but §2's table is not

Supplementary Table 5 lists first- and second-degree pairs. Among the candidates it gives exactly
one male–male first-degree pair: **`Aes12`–`Aes19`, `pi_HAT` 0.427, 272,655 SNPs, `same Y HG = yes`**
— precisely and only the merge registered in §8. The pooling rule needed no revision.

**But the wider family structure in §2's table is the compilation's, not the publication's.** Table 5
contains no entry pairing `Aes1`, `Aes20`, `Aes21` or `Aes23` with anything. The compilation
describes them as "2nd or 3rd degree" relatives within a twelve-member Family D; Table 5's stated
scope is first and second degree, so a genuinely third-degree pair would be absent legitimately —
but a second-degree one should have appeared. The compilation's kinship notes are therefore a
re-derivation of unstated provenance, exactly like its `Y-Haplotree-Variant` calls, and §2's "seven
kin groups" is **sourced to the compilation and must be attributed to it**, not to Furtwängler 2020.

The consequence runs in the cohort's favour and is deliberately not banked: on the publication's
evidence alone there are **14 observation units with no documented relationship among them**, rather
than seven kin groups. That would make Test C stronger than registered. It is not claimed, because
absence of a pair from a table scoped to first and second degree is not evidence of unrelatedness,
and because neither source measures the thing that actually matters — patrilineal sharing beyond
third degree, which is unmeasured in both. **FC3 stands unchanged and is still the right guard.**

Two further pairs are first-degree but do not reduce the count: `Aes11`–`Aes17` and `Aes14`–`Aes15`
each pair a candidate with a female (`Aes11`, `Aes15`), so no Y is shared and no merge follows.

### Two observations recorded in passing, neither bearing on this test

- **`Aes12`–`Aes24` are first-degree with `same Y HG = same clade`**, yet `Aes12` is published
  `PF3239` and `Aes24` published `Z6488`. Two men who must carry the same Y chromosome carry
  published labels one node apart. This is a fresh instance of the point already made from Table 5
  for Oberbipp: **the published terminal is a coverage floor, not a placement.** `Aes24` is not in
  the cohort and nothing here depends on it.
- **`SX11` is internally inconsistent in Table 1**, carrying terminal mutation `PF3239` against
  `YHG = G2a2b2a1a1`. `PF3239` is `G2a2a1a2a1`; the two cannot both be right. `SX11` is Niederried
  Ursisbalm and is the single non-Aesch/Oberbipp/Muttenz individual in the published `PF3239` tally
  of 20, so any future count over that tally should treat it as unresolved.

---

## Registered outcome, 2026-07-27

Applied mechanically from §§5, 6 and 9. Interpretation is deliberately not in this section.

### Callability against the registered prediction

| quantity | predicted (§4) | observed |
|---|---|---|
| individuals callable at `Z6219` | 4.0 | **5** |
| individuals callable at `L166` | 5.9 | **12** |
| individuals callable at **both** | 2.1 | **4** |

The model under-predicted throughout, in the direction the `NRY` defect implies: `Aesch13` alone was
registered at essentially zero power and is callable at both.

### The four individuals callable at both

| individual | `Z6219` | `L166` |
|---|---|---|
| `Aesch12` | DERIVED (2/0) | ancestral (0/2) |
| `Aesch13` | DERIVED (2/0) | ancestral (0/2) |
| `Aesch23` | DERIVED (2/0) | ancestral (0/6) |
| `SNPRA62` | DERIVED (2/0) | ancestral (0/2) |

`SNPRA61` is DERIVED at `Z6219` (2/0) with a single ancestral read at `L166`, so it is callable at
`Z6219` only. Under the registered pooling rule `Aesch12` and `Aesch19` are one unit; `Aesch19` has no
`Z6219` coverage and 5 ancestral reads at `L166`, so the pooled unit is unchanged in direction.

Across all 15 candidates: **zero `Z6219`-ancestral reads** and **zero `L166`-derived reads**.
`L166` drew 45 ancestral reads, `L167` 24 ancestral against 1 derived, `Z6219` 15 derived.

### Damage evidence, required by §8 for decisive C>T calls

`Z6219` is C>T at chrY:13782251. All ten derived reads carry BQ 32-41 and MAPQ 37. Eight of ten sit
17-51 bp from the 5' end. The two exceptions are both in `SNPRA62`, at 1 and 4 bp. The libraries are
UDG-treated and their measured chrY damage is low — terminal C>T of 1.0-1.5%, falling to ~0.15% at
position 1 and ~0.05% by position 4 — so even the worst-placed pair has a joint damage probability
below 1e-5. Damage is excluded for all five.

### Falsifiers and predictions

- **FC1 — not triggered.** It requires two or more independent kin groups callable at both and
  `Z6219`-ancestral. Zero `Z6219`-ancestral reads exist in the cohort.
- **FC2 — not triggered.** It requires an individual `Z6219`-ancestral and `L166`-derived. Neither
  allele appears anywhere: zero ancestral reads at `Z6219`, zero derived reads at `L166`.
- **FC3 — cannot be evaluated as written, and the rule is degenerate.** See below.
- **PC1 — met.** Kin groups outside Oberbipp are `Z6219`-derived and `L166`-ancestral.
- **PC2 — no split.** All five callable at `Z6219` are derived. §6 registered uniformity in either
  direction as **uninformative rather than confirmatory**, and it is recorded that way. This is not
  evidence for H1.
- **PC3 — holds.** The ten YFull positions that drew zero reads across 15 Oberbipp libraries drew
  zero reads here: `C101821`, `FGC5721`, `FT191098`, `FT91632`, `Z6206`, `Z6208`, `Z6220`, `Z6287`,
  `Z6499`, `Z6519`. The mapping control passes.

### FC3 is degenerate and this is a defect in the registration

FC3 fires if "every callable candidate turns out to carry an identical Y haplotype across the backbone
control set". No candidate differs from any other at any backbone marker: every informative cell is
derived and there are no ancestral calls. On a literal reading FC3 obtains and §9 sends the outcome
to H0.

That reading is not usable, and the reason is a flaw in the falsifier rather than a fact about the
cohort. The backbone control set **is** the published `G2a2a1a2a1` path, and the cohort was selected
on a label asserting membership of that clade. A set of markers defining a clade cannot discriminate
within it. FC3 therefore fires with probability 1 for any cohort this document could have assembled,
and a falsifier that cannot fail to fire tests nothing.

It is also true that identity is not established: 95 of the 150 backbone cells are `no_coverage`, so
what is observed is absence of contradiction, not demonstrated sameness.

**Recorded as: FC3 not evaluable.** The outcome is *not* being routed to H0 on it, and equally the
independence it was written to protect is *not* claimed. The registration is not amended to repair
FC3 — under `iceman-y/PROTOCOL_extending_analyses.md` a rule found wanting after seeing its effect gets its own
registration and a re-run, never a patch applied to the run that exposed it.

### Outcome under §9

**PC1 met → H1 supported at a second site**, with the number stated rather than the word
"replicated":

- **4 individuals** callable at both markers, all `Z6219`-derived and `L166`-ancestral;
- **3 documented kin groups** — Aesch Family D (`Aesch12`+`Aesch19`, `Aesch23`), Aesch Family B
  (`Aesch13`), and Muttenz (`SNPRA62`, plus `SNPRA61` at `Z6219` only);
- **2 sites** beyond Oberbipp.

**These are not demonstrably independent patrilines and this design cannot make them so.** §8's
independence rule stands: kin groups negative at third degree are not thereby Y-independent, and
neither Furtwängler 2020 nor the compilation measures patrilineal sharing beyond that. The honest
count of demonstrably independent Y lineages carrying this state remains **1**, now observed at three
sites rather than one.

### Post-hoc, labelled as such: `Z6135`

Not registered as a decisive marker and not predicted. `Z6135` (C>T, chrY:7657979) drew **zero reads
in all 15 Oberbipp libraries** and was `uninformative` in Test B; here it is covered in 11 of 15 with
**23 derived reads and no ancestral read**, 8 of them callable at depth. Read evidence on the three
deepest libraries puts the derived reads at 1-49 bp from the 5' end with BQ 41 and MAPQ 37.

Its behaviour matches `Z6219` — derived where `L166`, `L167` and `Z6134` are ancestral — in men who
are callable at `L166`, including `Aesch1`, `Aesch4`, `Aesch20`, `Aesch21` and `Aesch22`, none of
which is callable at `Z6219`. `Z6209` (G>A) is derived in `Aesch22` at depth, with the same
orientation.

This is a scan result over 22 positions and is reported as post-hoc, not as a prediction that came
true. What it justifies is a *new registration*, not an extension of this one: the question of whether
`Z6219`, `Z6135` and `Z6209` form a block below `PF3239` and above `L166` is a different question from
the one this document asked, and the markers were chosen after seeing which ones had coverage.
