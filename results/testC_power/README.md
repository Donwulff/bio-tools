# Test C power estimate — inputs and outputs

Produced 2026-07-26 for `PREREG_testC_aesch_muttenz.md`. Commands in `RUNLOG.md`.

**No Aesch or Muttenz sequence data exists in this repository or on this machine.** Everything here
was computed from published read counts, an external coverage proxy, and marker rates measured on
the Oberbipp cohort in `results/swiss15/` and `results/testB/`. These are the numbers the test was
registered against, committed before it ran so that the registered expectation can be checked
against the outcome rather than recalled.

## Files

| file | what it is |
|---|---|
| `candidates.tsv` | the 15 hedged `G-L166*` individuals at Aesch and Muttenz, with the `NRY` coverage proxy |
| `candidates_lineage.tsv` | the same cohort as **observation units**, after the registered pooling rule |
| `proxy_calibration.tsv` | `NRY` for the Oberbipp/Rapperswil samples, used to fit proxy → chrY reads |
| `testC_calibration.tsv` | fitted rates: proxy ratio per sample, and per-marker read rate with Poisson interval |
| `testC_candidates.tsv` | per-candidate λ and P(callable) at `Z6219`, `L166`, `FGC5671`, `L167` |
| `testC_power.tsv` | cohort-level expectations and split power, four markers |
| `testC_lineage_*.tsv` | the same, on the lineage-collapsed cohort — **these are the registered numbers** |

## Which table is the registered one

`testC_lineage_power.tsv`. The four-marker `testC_power.tsv` is kept because `RUNLOG.md` reproduces
it and because it carries the `FGC5671` and `L167` rates, but its final row — the joint over all
four markers, 0.05 expected — is arithmetic, not a plan. Nobody proposed requiring four
simultaneous calls. The quantity that matters is the **`Z6219`+`L166` pair**, because the
intermediate-branch signature is derived at one and ancestral at the other **in the same man**.

## Reading `candidates_lineage.tsv`

The rows are observation units, not individuals. Only `Aesch12+Aesch19` is a merge, because that is
the only documented 1st-degree male–male pair among the candidates and therefore the only pair
guaranteed to share a Y chromosome by descent. Every other row is a single individual even where
the `note` column records 2nd- or 3rd-degree kinship, because such a link may run through a female
and pooling across an unshared Y would manufacture a false split inside one unit.

The `note` column preserves the kin group so the units can be re-collapsed if the backbone control
later shows the families are Y-identical after all. That outcome is registered as **FC3** and its
consequence is H0, not support for anything.

## The number most likely to be misread

`expected_callable` for `Z6219` is **4.0 of 15**, and it is a point estimate from a rate fitted on
**ten reads**. Its Poisson interval puts the same quantity between **1.3** and **7.9**. The interval
is not a caveat attached to the estimate; at these counts it *is* the estimate. Any statement about
Test C's power that quotes 4.0 without the bracket is overstating what was known before the run.

## `Aesch13`

Carries `NRY = 17` and a predicted chrY yield of about 11 reads. It is registered as a candidate and
expected to return no power. It is listed rather than dropped so that reporting it as uninformative
is a recorded outcome instead of a silent exclusion.
