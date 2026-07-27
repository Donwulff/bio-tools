# Protocol: extending an analysis past its pre-registered tests

Written 2026-07-26, after Tests A and B of `iceman-y/prereg/Z6219_node.md` returned results that immediately
suggested three follow-up tests none of which were registered. That situation is the normal one, not
an exception, and it needs a rule rather than a judgement call each time.

## Why this exists, given that we compute no p-values

Haplotree work mostly resolves unambiguously or not at all, so the multiple-comparisons framing
imports badly. Three specific failure modes survive the translation, and they are what this protocol
guards:

1. **Marker shopping.** Scanning 22 positions and reporting the one that fits is the same error as
   scanning 22 hypotheses for p < 0.05, with no p-value anywhere in sight. Test B found exactly one
   co-segregating marker out of 22. Whether that is a finding or a coincidence depends entirely on
   whether the counting rule was fixed in advance.
2. **Silent conversion of "untested" into "negative".** 17 of Test B's 22 positions were
   uninformative. Nothing but a written rule stops that being reported as "no other marker
   co-segregates", which is a different and much stronger claim.
3. **Choosing the statistic after seeing the data.** At `FGC5671` the first look was "four of five
   reads sit near the 3' terminus" — which suggests damage. Weighting by the library's own damage
   rate reverses it. Both statistics were available; only one was specified in advance. The rule that
   saved this was §8 of the pre-registration requiring *both*, written before the data existed.

The benefit is therefore not statistical hygiene in the p-value sense. It is that **the criterion is
fixed while it can still go either way**, and the record afterwards distinguishes what was predicted
from what was noticed.

## The trigger: when does a test need registering?

One question, applied before running anything:

> **Can I name, right now, an outcome of this test that I would be disappointed by?**

- **Yes → register it first.** Disappointment is the signal that a preferred result exists, and a
  preferred result is what bends interpretation.
- **No → just run it and log it.** Instrument checks, tool validation, coverage counts, "does the
  index build" — these have no side to favour.

The question is deliberately about *wanting*, not about *importance*. A test can be trivial and still
need registering; a test can be central and not need it. `bwa index` completing is central to the
work and needs nothing. `FGC5671` in one more sample is a small addition and needs registering,
because there is an outcome I would rather have.

### Second trigger, for cases where the first is unclear

> **Would this test's result change something already written down as a conclusion?**

If yes, register. This catches tests that feel neutral but function as re-litigation of a settled
claim — the direction from which F6 arrived.

## The lightweight form

A registration is not a document. For most extensions it is five lines appended to the relevant
prereg under `## Extensions`, and it must contain all five:

```
### E<n> — <one-line question>            [date]
Prior state:  what data bearing on this I have already seen, named explicitly.
Prediction:   what I expect, and what the opposite outcome would mean.
No-power:     what result counts as "could not tell", stated before running.
Decision:     what changes in the repository depending on each outcome.
Status:       REGISTERED / RUN <date> / ABANDONED <reason>
```

`Prior state` is the peeking declaration and is the line most worth being pedantic about. "I have
read the Oberbipp and Sardinian calls at this marker but not the Iceman's" is the sentence that makes
everything after it interpretable.

`No-power` is not optional and is never merged into either outcome. It is the most likely result at
ancient-DNA depths and must be reportable as itself.

## Amendment, extension, or new registration

| situation | what to write |
|---|---|
| A registered test's rules turn out to be wrong or incomplete | **Amendment** to the existing prereg, stating what it threatens |
| A new test of an existing hypothesis, new sample or marker | **Extension** (`E<n>`), five-line form |
| A new hypothesis, or a cohort with its own power question | **New pre-registration** document |

The distinction that matters is the middle row: extensions are cheap on purpose, so that the excuse
"registering this is more work than running it" never becomes true.

## Rules that apply to all three

**The ratchet.** Once registered, a test's result is recorded whatever it says. A registered test
that was not run stays on the outstanding list until it is run or explicitly abandoned with a reason.
`RUNLOG.md`'s "Not done, and explicitly outstanding" section is where that lives.

**Post-hoc is a label, not a verdict.** Results found by noticing rather than predicting are worth
reporting; they are simply reported as post-hoc. Suppressing them is worse than labelling them. What
is prohibited is presenting one as though it had been predicted.

**No rule change after seeing its effect.** If a filter, threshold or call rule turns out to be
wrong, the repair gets its own registration and a re-run of every affected dataset — it is never
applied to the run that revealed it. This is why `MAPQ >= 25` is still in force despite being known
allele-biased, and why `uniqueness_audit()` reports without acting.

**Counting rules are fixed before the scan.** For any scan over markers or samples, write down
before running what counts as tested, untested, positive and negative. Test B's
`uninformative` / `stays_with_block` / `splits_with_X` split is the worked example.

## When to skip all of this

Registration is for tests of claims. It is not for:

- building indexes, staging data, checking file provenance;
- tool validation against a known answer (the `(H)`-flag reproduction);
- re-running a committed pipeline to check byte-identity;
- reading a table to find out whether a test is even possible, *provided* the columns bearing on the
  hypothesis are not read — and if that cannot be separated, register first.

The last one has teeth: checking whether `I5118` has coverage at `FGC5671` is fine; reading which
allele he carries is the test.

## Current queue

| id | test | why it needs registering | status |
|---|---|---|---|
| E1 | `FGC5671` in the Iceman | Would resolve whether the second `Z6219`-block SNP is real or a private Oberbipp variant. I want it derived. | **RUN 2026-07-26**, prediction held; F3 does not obtain |
| E2 | `CGG017683` (Crimea), genotyped at the YFull `L166` positions | Asymmetric falsifier: could exclude the sample, could not support it. | **RUN 2026-07-26**, H0 — zero reads at every decisive position; closed as depth-limited |
| E3 | Uniqueness filtering (`XT:A:U`, `X0=1`) replacing `MAPQ >= 25` | Changes a rule governing every call; promoted to its own document, `iceman-y/prereg/uniqueness_filter.md` | **RUN 2026-07-26, NOT ADOPTED** — U1 and U2 both fail; it rescues the paralogue at `FGC5687` |
| C | `Z6219` and `L166` in the Aesch/Muttenz cohort | A cohort with its own power question, so a new document rather than an extension: `iceman-y/prereg/testC_aesch_muttenz.md` | **RUN 2026-07-27. PC1 met**, FC1/FC2 not triggered, PC3 holds; FC3 not evaluable and left unpatched |
| D | Do `Z6219`, `Z6135` and `Z6209` form a block below `PF3239` and above `L166`? | The markers were chosen after seeing which had coverage in Test C, so this is a new hypothesis and not an extension of it. Needs its own power statement and its own falsifiers | **NOT REGISTERED.** Post-hoc observation only |
