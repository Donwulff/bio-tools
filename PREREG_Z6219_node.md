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

---

## Amendment 1 (2026-07-26): allele-aware mappability, registered before running

**Prompted by two external signals**, reported by the user and treated here as data rather than as
verified fact: **FTDNA does not carry `Z6219` at all**, and **YFull rates it 3 of 5**, among the
least reliable in its set. The precise meaning of YFull's scale has not been verified here and is
not relied on. What matters is that two independent curators have either excluded or down-rated this
marker, which is a reason to look for an artifact this project has not yet tested for.

**A limitation of the evidence in §2 is acknowledged.** The mappability sweep
(`results/mappability/`) tiled reads cut **from the reference**, so every read carried the
**ancestral** allele. `Z6219 frac_recovered = 1.000` therefore establishes only that an *ancestral*
read returns. It says nothing about a read carrying the derived `T`, which has a lower alignment
score against this locus and may score better against a paralogous one. §2's claim that "mappability
is not a confound" is hereby narrowed to ancestral reads and must not be quoted more broadly.

### F5 — allele asymmetry / paralogy (new falsifier)

`annotate/y_mappability.py` will be extended with an `--allele` option that substitutes the derived
base at the marker position before tiling, and run for both alleles at `Z6219` and, as controls, at
`PF3239`, `L166` and `Z6287`.

**F5 obtains if**, for `Z6219`, derived-allele reads recover materially worse than ancestral-allele
reads, **or** if any derived-allele read maps to a locus other than `chrY:13782251`. Either shows a
near-identical paralog capable of donating false derived reads, which is a complete alternative
explanation for the Oberbipp observation and would falsify H1 as currently supported.

**Threshold, fixed now:** a drop in `frac_recovered` of **> 0.10** between alleles at 45 bp, or
**any** off-target placement of a derived-allele read at MAPQ >= 25, triggers F5. The same
comparison at `PF3239`, `L166` and `Z6287` is the control: if all four markers show similar
asymmetry, it is a general property of short-read alignment and not specific to `Z6219`.

**Prediction if H1 survives:** `Z6219` behaves like the control markers, and the curators'
reservations are explained by something other than paralogy at this position — for example depth or
region annotation — which would then need naming rather than assuming.

**What F5 does not cover.** A paralog that is *absent from the reference entirely* cannot be
detected this way, since reads from it have nowhere correct to go. That failure mode remains
untested and must be stated as a residual risk regardless of outcome.

**Not yet read at time of writing:** any allele-substituted mappability result, and the per-sample
depth profile around `chrY:13782251`.

## Amendment 2 (2026-07-26): MAPQ threshold is allele-biased — registered before checking

F5 was diagnosed to its mechanism and it is **not** paralogy. Single-read test against the no-alt
reference, `bwa aln -n 0.01 -k 2 -l 1024`:

    L166_anc    MAPQ=37  XT=U X0=1 X1=0 NM=0   PASS
    L166_der    MAPQ=23  XT=U X0=1 X1=1 NM=1   FILTERED  <-- unique, and discarded
    Z6208_der   MAPQ=20  XT=U X0=1 X1=2 NM=1   FILTERED  <-- unique, and discarded
    Z6219_der   MAPQ=37  XT=U X0=1 X1=0 NM=1   PASS
    PF3239_der  MAPQ=37  XT=U X0=1 X1=0 NM=1   PASS

Every alignment is unique (`XT:A:U`, `X0=1`). MAPQ drops only because a *suboptimal* hit exists
elsewhere and the read carries its one mismatch. **Our `MAPQ >= 25` threshold therefore discards
uniquely-mapped derived reads at `L166` and `Z6208`, and keeps every ancestral read.** The bias is
allele-specific and its direction is toward calling ancestral. Identical in `working`, `noalt` and
`hs37d5`, so it is a property of the locus and the aligner, not of the custom reference.

### F6 — the Oberbipp ancestral calls may be incomplete

**This threatens a conclusion already published in this repository.** Oberbipp is recorded as
`L166` ancestral 11/11, 0 derived, which is the basis for refuting the circulating "7 of 10" claim.
If those libraries contain derived `L166` reads sitting at MAPQ 20-24, the true call is **mixed**,
not ancestral, and the refutation must be restated.

**Check, registered now:** re-examine every read at `chrY:21843737` in all 15 Swiss BAMs with **no
MAPQ floor**, reporting allele, MAPQ, `XT`, `X0`, `X1` and read-terminus position for each.

- **F6 obtains** if any Oberbipp sample carries one or more `L166`-derived reads with `XT:A:U`,
  `X0=1` and MAPQ >= 20 that the 25 threshold excluded. The ancestral calls are then incomplete and
  the "7 of 10" refutation is reduced accordingly.
- **F6 does not obtain** if the only sub-threshold reads at that position are `XT:A:R`/`X0>1`
  (genuinely ambiguous) or absent. The ancestral calls then stand as recorded.

**The threshold will not be changed as a consequence of this result.** Switching from `MAPQ >= 25`
to a uniqueness criterion (`XT:A:U` and `X0=1`) is the principled repair and is very likely correct,
but adopting it *after* discovering that it alters a specific finding is exactly the post-hoc rescue
this project's protocol exists to prevent. It is recorded here as a proposal requiring its own
registration and a re-run of every affected dataset, not as a change to be made now.

**Consequences that follow regardless of F6**, since they are arithmetic rather than findings:
derived-read recovery at `L166` is 0.156 at 45 bp and at `Z6208` 0.044, so a genuinely derived
individual loses most of his reads at those sites. `I5118`'s single derived `L166` read implies of
the order of six underlying molecules, and any sample reported "no coverage at `L166`" may be a
derived man whose reads were filtered. Every `L166`/`Z6208` derived count in this repository is a
**lower bound**, and no `no_coverage` at those two sites may be read as evidence of anything.

**Not yet read at time of writing:** the sub-threshold read content at `chrY:21843737` in any sample.

## Amendment 3 (2026-07-26): corrections to §7 power and to the Amendment 1 premise

**The `Z6219` rating was understated, not overstated.** Amendment 1 recorded YFull's 3/5 and I then
noted that YFull's FAQ deems 2–5 stars "good quality" and uses all of them to build the YTree,
treating the signal as unremarkable. That reply answered the absolute scale and missed the point
actually being made: **`Z6219` ranks 26th of the 32 SNPs YFull uses to define `L166`**, ordered by
reliability, in a set that includes two 1-star markers. Bottom-quartile within its own block is a
meaningful position regardless of where 3 stars sits on the global scale. The reservation was real.

It does not change the F5 outcome — `Z6219` is 1.000 at both alleles at 35, 45, 60, 100 and 150 bp,
with zero MQ0 and zero off-target placements, in all references tested. Whatever YFull's rank
reflects, this project has not reproduced it, and that gap is now itself unexplained rather than
dismissed.

**§7 Test B power is revised upward.** That section estimated ~20 distinct catalogued positions at
the `G-L166` node, of which 9 are already tested, leaving ~11. **YFull defines `L166` with 32 SNPs**,
so our catalogue under-represents the block and the co-segregation scan has more to work with than
registered. Revised before the scan is run: the marker set for Test B should be sourced from YFull's
`L166` definition rather than from `resources/marker_index.tsv.gz` alone, and the F3 threshold
("no other marker co-segregates") must be evaluated over that larger set. Obtaining that list is a
prerequisite to Test B and has not been done.

No result from Test A, Test B or the F6 recheck influenced this amendment; F6 was the only one run
at time of writing and it returned "does not obtain".

---

## Extensions

Registered under `PROTOCOL_extending_analyses.md`. Each is a test of a hypothesis already stated
above, using the five-line form; the rules of §8 govern all of them unchanged.

### E1 — Is `FGC5671` derived in the Iceman? [2026-07-26]

**Prior state.** I have read every `FGC5671` value this project has produced: the 15 Swiss samples
in `results/testB/swiss_yfull_L166_defining.tsv` (`MX210` DERIVED 2 derived / 0 ancestral; `MX213`,
`SX10`, `MX209` each `nocall_damage_prone_1read`, derived direction; the remaining 11 `no_coverage`),
both Sardinians in `results/testB_unhedged/` (`I14677` 1 ancestral, `I14678` 3 ancestral), the Test B
verdict `splits_with_Z6219`, the allele-aware mappability (1.000 at both alleles at 35/45/60 bp,
zero MQ0), and the read-terminus and damage-profile forensics at the site. I have **not** looked at
any Iceman read, pileup, call or coverage figure at `chrY:7,784,648`. The marker is absent from
`markers/L166_defining.txt`, so it was not in the 9-marker set of the driver validation, and no
committed Iceman table contains the position or the name — checked by grep before writing this.
The catalogue holds exactly one name at that coordinate, so no synonym could have leaked it.

**Prediction.** Derived. The Iceman is derived at `PF3239`, `Z6219` and `L166` alike (§10), so a
marker lying anywhere on the shared path down to the `L166` block must be derived in him. This is
the outcome I want, which is why the test is registered rather than simply run.

- **Derived** → `FGC5671` is a shared branch marker, not a private Oberbipp variant. **F3 does not
  obtain**: a second position co-segregates with `Z6219`, and homoplasy (H3) would have to have
  struck twice in the same lineages, at two independent positions, to reproduce it.
- **Ancestral (>= 2 reads)** → an `L166`-derived man lacks `FGC5671`, so it cannot be an `L166`
  block marker and **YFull's block assignment for it is wrong**. Test B's co-segregation count
  falls to one (`Z6219` alone), **F3 obtains**, and the node is reported as unresolved between H1
  and H3. Ancestral reads at this G>A site are damage-robust by direction, so two suffice, by the
  same argument as F1.
- Either way this says nothing about the Iceman's own placement, per §10.

**No-power.** Any of `no_coverage`, `nocall_noreads`, `nocall_damage_prone_1read`,
`low_power_1read_ancestral`, or a mixed call failing `site_qc` / `pct_mq0 < 30%`. Reported as H0
for this extension. Given that 11 of 15 Oberbipp samples have no coverage here, H0 is a live
outcome and not a fallback.

**Decision.** `NOTES.md` gains the result and `RUNLOG.md` the command, whichever way it lands. The
"the gap that decides it" paragraph in `NOTES.md` is rewritten to state what was found, not
removed. If ancestral, the Test B section is amended in place to record that its
`splits_with_Z6219 = 2` overcounted, and the F3 verdict changes.

**Fixed in advance — one BAM.** The test runs against
`iceman.oetzi.UDG_merge_combined.mapped_rmdup.pair.prim_rmdup.sort_rmdup.coord.bam`, the merged BAM
every other Iceman finding in this repository rests on, and nothing else. There are ~20 further
Iceman BAMs on disk from earlier assemblies and reference builds. If the canonical BAM gives no
coverage, that is H0; trying another library until one yields a read is exactly the search this
protocol exists to prevent, and requires its own registration.

**Secondary scan, registered because it cannot be avoided.** The tool emits all 22 testable
positions of `markers/yfull_L166_defining.txt` in one pass, so I will see the other 21 whether or
not I intend to. Their meaning is therefore fixed here, before the run:

- 9 of the 22 are `markers/L166_defining.txt` and are already known derived in the Iceman (RUNLOG,
  driver validation). They are **excluded from the count** as non-independent and serve only as a
  positive control: anything other than 9/9 derived invalidates the run.
- Of the remaining 13, a position is **block-confirmed** if the Iceman is `DERIVED`,
  **block-refuted** if he is `ancestral` under the registered rules, and **untested** otherwise.
  `untested` is never counted as either.
- Prediction: mostly block-confirmed, and `block-refuted > 0` is the interesting outcome — it would
  mean YFull's 32-SNP `L166` definition contains positions that an `L166`-derived man does not
  carry. That claim, if it arises, is registered here and is not post-hoc.
- No marker outside these 22 is tested under E1.

**Status.** REGISTERED 2026-07-26.
