# Pre-registration: should `MAPQ >= 25` be replaced by a uniqueness criterion?

Written 2026-07-26, **before** the measurement below was run. This is E3 of
`iceman-y/PROTOCOL_extending_analyses.md`, promoted from an extension to its own document because it is not
a test of any hypothesis in `iceman-y/prereg/Z6219_node.md` — it is a change to a rule that governs **every**
call this repository has ever made, and it therefore needs its own scope and its own decision rule.

Amendment 2 of `iceman-y/prereg/Z6219_node.md` registered this as a proposal and explicitly refused to act on
it: *"adopting it after discovering that it alters a specific finding is exactly the post-hoc rescue
this project's protocol exists to prevent."* This document exists to make the adoption decision
**before** the effect on any finding is known.

---

## 1. The defect being repaired

`bwa aln` MAPQ conflates two different things: whether a read was placed at exactly one locus, and
whether it matched perfectly there. Measured against the no-alt reference at `-n 0.01 -k 2 -l 1024`:

    L166_anc    MAPQ=37  XT=U X0=1 X1=0 NM=0   PASS
    L166_der    MAPQ=23  XT=U X0=1 X1=1 NM=1   FILTERED  <-- unique, and discarded
    Z6208_der   MAPQ=20  XT=U X0=1 X1=2 NM=1   FILTERED  <-- unique, and discarded
    Z6219_der   MAPQ=37  XT=U X0=1 X1=0 NM=1   PASS

Every one of these alignments is unique. MAPQ falls only because a *suboptimal* hit exists elsewhere
and the read carries its one mismatch. A derived read carries a mismatch **by definition**, so
wherever `X1 >= 1` the threshold discards derived reads and keeps ancestral ones. The bias is
allele-specific, its direction is toward calling ancestral, and it is identical in the `working`,
`noalt` and `hs37d5` references, so it is a property of the locus and the aligner.

**The proposed repair:** accept a read iff `XT:A:U` **and** `X0 == 1`, regardless of MAPQ.

## 2. Why the repair is not obviously safe

This must be stated before the measurement, because it is the reason a decision rule is needed at
all rather than just a patch.

`MAPQ >= 25` is not only an allele filter. It is also the only thing standing between the caller and
reads that have a plausible alternative home. Admitting every `X0 == 1` read admits reads with
`X1 >= 1` — reads for which the aligner *found* another candidate locus and preferred this one by a
margin the new criterion does not examine. At a paralogous or X-homologous position a damaged read
can be uniquely *best* in the wrong place.

So the repair has its own failure mode, and it runs in the **opposite direction — toward calling
derived.** The central claim of this entire project is a derived call. A filter change that quietly
manufactures derived reads would produce exactly the result I want, by construction, and would be
almost impossible to detect afterwards. **U2 below exists solely to catch this**, and it is a hard
fail with no override.

## 3. Hypotheses

**Hu — the repair works.** `XT:A:U && X0 == 1` is allele-neutral where `MAPQ >= 25` is not, and
admits no reads that belong elsewhere.

**Hb — the repair trades one bias for another.** It removes the ancestral bias and introduces a
derived one at exactly the positions where paralogy is real.

**H0 — no power.** The measurement cannot distinguish them.

## 4. The measurement, and why it is not circular

Everything the decision rests on is measured on **simulated reads cut from the reference itself** —
exact copies, no sequencing error, no damage, no biology. `y_mappability.py` already works this way:
for each marker it cuts every read of length L overlapping the position and maps it back with the
pipeline's own aligner and options, once with the ancestral base and once with the derived base
substituted.

This is what makes the test non-circular. **The truth is known by construction.** Both allele tiles
are exact reference sequence, so *any* difference in recovery between them is produced by the filter
and nothing else. A correct filter must return the same number for both. No sample's genotype, and
no finding of this project, enters the measurement.

**Marker set, fixed now:** `markers/L166_defining.txt` (9), `markers/Z6494_exclusion.txt` (3) and
`markers/backbone_control.txt` (10) — the three established sets. They are chosen because they
contain, in advance and for independent reasons: the two markers known to be asymmetric (`L166`,
`Z6208`); the marker known to be unrecoverable because its reads genuinely belong elsewhere
(`FGC5687`, 0.000 recovery with 45/45 MQ0 and 29 off-target placements); and the ten backbone
positives that every established result depends on. Read lengths 35, 45, 60, 100. Targets: the
`working` reference and `noalt`.

**Tooling change:** `y_mappability.py` gains `n_unique` and `frac_recovered_unique` **as additional
reported columns**. The existing `n_mq_ge_min` and `frac_recovered` columns are computed exactly as
before and their values do not move. The new filter is measured, not applied — the same
report-without-acting discipline `ylib.uniqueness_audit()` already follows. No call in the
repository changes as a result of running this measurement.

## 5. Adoption criterion, fixed before the measurement

Uniqueness filtering is adopted **iff all three hold**. Any single failure means it is not adopted.

**U1 — the bias is actually repaired.** For every marker in the set, at every read length, the
absolute difference between ancestral and derived `frac_recovered_unique` must be **<= 0.05**, and
no marker's derived recovery may be *lower* under the uniqueness criterion than under MAPQ. The
threshold is 0.05 rather than 0 to absorb tile-edge effects at the ends of the simulated window; it
is not a tolerance for real asymmetry, since on exact reference tiles there is no such thing.

**U2 — the cure admits nothing it should not.** At every marker with `n_offtarget > 0` under the
current measurement — `FGC5687` above all — `frac_recovered_unique` must **not exceed**
`frac_recovered`. Any increase means the criterion is admitting reads the whole-genome measurement
says belong elsewhere. **Hard fail, no override**, for the reason given in §2.

**U3 — known answers survive.** Recomputed under the new criterion on the real BAMs, all three
established positive controls must reproduce exactly: Iceman **9/9 derived** at `L166_defining`,
**10/10 derived** at `backbone_control`, **3/3 ancestral** at `Z6494_exclusion`.

### The firewall

**The decision is made on U1, U2 and U3 alone.** Every one of them is measured either on simulated
reads whose truth is known by construction, or against a published answer fixed long before this
document.

**Whether the Oberbipp `L166` calls move, whether `I5118` gains reads, and whether any finding of
this project is strengthened or weakened, are measured and reported but do not enter the decision.**
If they did, this would be the post-hoc rescue Amendment 2 forbids, merely with more ceremony.

**No partial adoption.** The filter is not adopted "at `L166` but not elsewhere", nor "for
transversions only". Marker-level rule selection is rule shopping and is prohibited outright.

## 6. Peeking declaration

Already read, and quoted above because it motivated this document: the five single-read alignment
records in §1; that `L166` derived recovery is 0.156 and ancestral 0.956 at 45 bp; that `Z6208`
derived recovery is 0.044; that `FGC5687` is 0.000 with 45/45 MQ0 and 29 off-target; that `Z6215`
is 0.778 and `Z6519` 0.000.

**F6 has already been run and does not obtain.** All 14 reads at `chrY:21843737` across the 15 Swiss
BAMs at *no* MAPQ floor are ancestral `C`, MAPQ 37, `NM=0`, with zero derived reads at any MAPQ.
This is declared here because it substantially reduces the conflict of interest: **the Oberbipp
`L166` result cannot change under any filter**, since there are no sub-threshold reads there to
admit. What can still move is `I5118`'s derived count, the `Z6208` calls, and any position currently
reported `no_coverage` at `L166`/`Z6208`.

**Not yet read at time of writing:** any `XT`/`X0` value from the simulated tiles under either
allele; any recomputed call under the uniqueness criterion; the effect on any sample's genotype.

## 7. Decision rules

- **All three hold →** adopt. `ylib` switches to the uniqueness criterion, every affected dataset is
  re-run, and both the old and new tables are kept alongside a diff table naming every call that
  moved. The change is applied to *all* datasets in the same commit or to none.
- **Any fail →** not adopted. `MAPQ >= 25` stands, `uniqueness_audit()` continues to report without
  acting, and Amendment 2's consequence remains in force: every `L166`/`Z6208` derived count in this
  repository is a **lower bound**, and no `no_coverage` at those two positions is evidence of
  anything.
- **H0 →** if the measurement cannot separate the criteria — for instance if every marker recovers
  1.000 under both, leaving nothing to repair — report H0 and keep `MAPQ >= 25`. The absence of a
  demonstrated repair is not a reason to adopt an untested one.

## 8. What this will not claim

- Nothing about which filter is correct in general, for other aligners, or at other loci. The
  measurement covers 22 markers on chrY under one aligner with one option set.
- Nothing about reads the depositors already discarded. Those are gone from every dataset here and
  no filter of ours can recover them.
- If adopted, no re-interpretation of any finding is authorised in the same commit as the re-run.
  Changed calls are recorded first; what they mean is a separate question asked afterwards.

---

## Outcome (2026-07-26): NOT ADOPTED

Run the same day, after this document was committed. **U1 fails, U2 fails.** `MAPQ >= 25` stands
unchanged and `uniqueness_audit()` continues to report without acting.

- **U1 — fails.** 28 of 176 cells are asymmetric beyond 0.05 under MAPQ; **18** still are under the
  uniqueness criterion. No marker's derived recovery got worse, so the repair points the right way,
  but `L166` at 45 bp only moves 0.156 → 0.689 against a 0.05 threshold. Reduced, not repaired.
- **U2 — fails, on the control named in §5.** `FGC5687` goes 0.000 → **0.250** at 60 bp and
  0.160 → **0.890** at 100 bp, at a position with 23 and 5–6 off-target tile reads. Ancestral and
  derived move together, so it is not an allele effect: it is admission of reads whose alternative
  locus is real. 29 cells fail. Hard fail as registered, no override taken.
- **U3 — not run.** Adoption required all three; two had already failed.
- **Effect on this project's calls — deliberately not computed**, per the firewall in §5. See
  `NOTES.md` for why that is the safer choice rather than an omission.

The decisive observation, which was not predicted: uniqueness repairs the asymmetry **only at
>= 60 bp**, and admits the paralogue **only at >= 60 bp**. There is no read length at which it is
both safe and useful. Recorded as post-hoc.

§7's "any fail" branch therefore applies in full: every `L166`/`Z6208` derived count in this
repository remains a lower bound, and that caveat is now permanent rather than pending a repair.
