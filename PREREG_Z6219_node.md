# Pre-registration: is `Z6219` a real node between `PF3239` and `L166`?

Written 2026-07-26, **before** any of the tests below were run. Inherits by reference from
`PREREG_swiss_neolithic_L166.md`: the per-SNP output rule, the MAPQ 30% rejection threshold, the
no-go regions, the coverage-denominator rule, the pooling prohibition, and the rule that a
documented patriline counts as one observation. Where this document is silent, that one governs.

---

## 1. What is being tested

Every published tree places `Z6219` **inside the `L166` equivalence block**: our catalogue gives it
`isogg G2a2a1a2a1a` / `yfull G-L166`, the same node as `L166` itself. The hypothesis under test is
that this is wrong, and that the true order is

    PF3239  ->  Z6219  ->  L166

with `Z6219` marking a real node strictly between them.

## 2. Status of the evidence so far — post-hoc, and labelled as such

This hypothesis was **not** derived from a registered prediction. It emerged from data collected to
answer a different question, and the pre-registration for that dataset explicitly stated these
samples could not settle it. That claim was wrong because it assumed the samples' published labels
were correct — the very thing under test. The evidence to date, all post-hoc:

- **Above `L166`.** Oberbipp: 10 reads derived at `Z6219` across four independent lineages
  (`MX210`, `MX213`, `SX10`, and Family A pooled as one observation), while ancestral 11/11 at
  `L166`. `L166` is a transversion and therefore damage-immune in both directions. All 10 `Z6219`
  reads sit 14–72 bp from the 5' terminus with `NM=1`, in libraries measured at 1–2% terminal C>T
  falling to background by position 1 — so the derived calls are not deamination
  (`NOTES.md`, "Z6219 Read-Terminus Check").
- **Below `PF3239`.** Sardinians: `I14677` derived 5/5 at `PF3239` and ancestral 4/4 at `Z6219`;
  `I14678` 3/3 and 2/2. An ancestral call at a C>T site is damage-robust by direction, since
  deamination cannot manufacture the ancestral `C`.
- **Mappability is not a confound.** `Z6219` and `PF3239` both return `frac_recovered = 1.000` at
  45 bp in all four references tested (`results/mappability/`), and all reads at both sites are
  0% MQ0.

Strength of evidence is not a substitute for having said in advance what would count. Hence this
document, and hence its orientation toward **falsification** rather than confirmation.

## 3. Why it matters, and the limit of that

If `Z6219` is a real node, then an ancient individual who is `Z6219`-derived and `L166`-ancestral
will be reported as `L166` by **any tool built on a modern haplotree**, because `Z6219` sits inside
`L166`'s equivalence block in ISOGG, YFull and FTDNA alike. Such a lineage may well have no living
male-line descendants; equivalence blocks can only be split by a sample on the intermediate branch,
and modern trees have no way to find one that did not survive. That makes this a systematic
confounder for ancient `G-L166` calls generally, not a curiosity about one SNP.

**The limit, stated up front so it is not overclaimed later.** This would *not* explain the
mislabelled Sardinians documented in `NOTES.md`: `I14677` and `I14678` are ancestral at `Z6219` as
well as at `L166`, so their over-deep labels are a separate failure. At most this mechanism explains
a subset of ancient `G-L166` calls — those that are genuinely `Z6219`-derived — and the analysis
must not be written up as though it accounted for all of them.

## 4. Hypotheses, including the ones that produce the same pattern

**H1 — real node.** `Z6219` lies strictly between `PF3239` and `L166`.

**H2 — equivalence (the catalogue's position).** `Z6219` is `L166`-equivalent, and the observed
Oberbipp pattern is an artifact of mismapping, coordinate error, or contamination.

**H3 — recurrence / homoplasy.** The `Z6219` derived allele arose independently in the Oberbipp
lineage and does not mark shared descent. **This is the dangerous alternative**: it predicts exactly
the same read-level pattern as H1 and is not distinguished by any amount of additional depth at
`Z6219` itself. The prior for it is not negligible — `L166` is catalogued as
`G-L166 & J-Y29712` and `L167` as `G-L166 & I-Y92994 & A00`, so **documented homoplasy already
exists at the neighbouring positions of this very block.**

**H0 — no power.** The tests below return no informative call. A permitted and expected outcome,
to be reported as H0 and never as support for either side.

## 5. Falsifiers, registered first

**F1 (primary).** Any individual **derived at `L166` and ancestral at `Z6219`**. Under H1,
`L166`-derived entails `Z6219`-derived, so a single well-supported counterexample falsifies the
ordering. Threshold: a derived `L166` call under the inherited rules, together with **>= 2** reads
ancestral at `Z6219` (two suffice because an ancestral `C` at a C>T site cannot be manufactured by
deamination), both sites at `site_qc pass` and `pct_mq0 < 30%`.

**F2.** `Z6219`-derived individuals who are **ancestral at `PF3239`**, or who fall outside
haplogroup G. Either indicates recurrence (H3) rather than a node on this branch.

**F3 (downgrade, not falsification).** If an exhaustive scan of catalogued positions at the
`G-L166` node finds **no** other marker co-segregating with `Z6219`, the node rests on a single SNP.
That is precisely what recurrence looks like, and H1 must then be reported as unresolved between H1
and H3 rather than as supported.

**F4.** If the 10 Oberbipp `Z6219` reads are not independent molecules — identical fragment
boundaries surviving deduplication — the observation collapses to fewer than the stated four
lineages and the section above must be withdrawn.

## 6. Predictions if H1 is true

**P1.** `I5118` is derived at `L166` (already established). Under H1 it must be `Z6219`-derived or
uncalled; ancestral would trigger F1.

**P2.** Of the catalogued positions at the `G-L166` node, at least one other should show the
`Z6219` pattern — derived in Oberbipp, ancestral in `I14677`/`I14678` — rather than the `L166`
pattern.

**P3.** Published-`PF3239` males from Aesch and Muttenz should **split**, some `Z6219`-derived and
some ancestral. Uniformity in either direction is uninformative rather than confirmatory.

## 7. Data and power, fixed before looking

**Test A — the F1 check.** `Z6219` is already in `markers/L166_defining.txt`, so its call for every
sample in `results/unhedged/` and `results/swiss15/` **already exists in committed tables**. See the
declaration in §9.

**Test B — co-segregation scan.** The catalogue holds **25 entries at `yfull_node = G-L166`**,
approximately 20 distinct positions after collapsing synonyms (`Z6213`/`S19530`,
`Z6516`/`FGC5675`, `Z6499`/`FGC5712`, `Z1370`/`Z1370.2`). Of these, **8 are transversions** and
therefore damage-immune, 12 are deamination-prone, 5 are transitions in the safe direction.
`markers/L166_defining.txt` currently tests **9**. The scan therefore has roughly **11 untested
positions**, and its power is bounded by their coverage in the Oberbipp and Sardinian libraries,
which is unknown at the time of writing and may well be the binding constraint. **A scan that
returns no covered position is H0, not F3** — F3 requires positions to have been *tested and found
not to co-segregate*, not merely to have been uncovered.

**Test C — new samples.** Out of scope for this registration. The 21 hedged Aesch/Muttenz
individuals are named in P3 as the natural cohort, but selecting and staging them is a separate
exercise requiring its own power statement, and nothing in this document authorises it.

## 8. Decision rules

- Per-SNP output only. No terminal label is assigned to any sample by this analysis.
- The pre-existing rules govern all calls. Where a C>T or G>A site is decisive, the **read-terminus
  and library damage-profile evidence** (`annotate/y_read_evidence.py`) must be reported alongside
  the count — the registered call rules judge damage by site class alone and have a known blind
  spot here, documented 2026-07-26. This is now protocol, not an ad-hoc rescue applied when a
  result depends on it.
- Family A counts as one observation. `I14677`/`I14678` count as one observation.
- Outcome to be recorded as exactly one of: **F1/F2 falsified**, **F3 unresolved (H1 vs H3)**,
  **H1 supported**, or **H0 no power**. "H1 supported" is not available if F3 obtains.

## 9. Peeking declaration

The `Z6219` genotype for `I5118`, and for every other sample in `results/unhedged/`, is present in
`results/unhedged/unhedged_L166_defining.tsv`, which is already committed. **I have not read those
values, and P1 is registered above without knowledge of them.** The Oberbipp and Sardinian `Z6219`
values *have* been read and are quoted in §2; they are the post-hoc observation that motivated this
document and are not evidence for it.

Test B has not been run in any form. The catalogue counts in §7 were obtained from
`resources/marker_index.tsv.gz` by mutation class and node only, with no sample data consulted.

## 10. What this analysis will not claim

- Nothing about geography, migration, or population replacement.
- Nothing about the Iceman's own terminal placement, which rests on separate evidence and is
  unaffected either way: he is derived at `PF3239`, `Z6219` and `L166` alike.
- Nothing about whether any lineage is extinct. Absence from modern trees is the neutral
  expectation for an arbitrary Neolithic patriline, not evidence of anything.
