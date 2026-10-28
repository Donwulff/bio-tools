# Ötzi's Y-haplogroup

Repository/workstream status, the unmerged Work assessments and interpretation limits are tracked in
[the dated situation review](../SITUATION.md). This page reports the project's read-level case;
metadata-only candidate research does not change its placements.

**Read-level evidence places the Iceman at `G-L166` under this project's local tree. `Z6208` is
derived, but is not treated as a separately resolved downstream node.**

Nature Communications 2025 (`s41467-025-61601-8`) reports `G2a2a1a2a1a1b (G-Z6208*)`. `Z6208` is a
real SNP and the Iceman is genuinely derived at it (DP 9). The catalogue/tree annotations inspected
for this analysis group it *inside* the `L166` block (`yfull_node=G-L166`), and the inspected ISOGG
longhand carries the provisional tilde (`G2a2a1a2a1a1b~`). The 30 Oberbipp, Aesch and Muttenz
libraries provide zero reads at Z6208, so they do not order it against L166. The local tree marks
the proposed separate node `refuted`; that is a project interpretation, not a refutation of the
derived allele. It also does not prove a universally basal terminal-star placement.

The local tree encodes two provisional splits with different evidential strength:

- **`G-Z6219` sits between `PF3239` and `L166`**, splitting the YFull `L166` block upward. Men from three
  Swiss Neolithic sites stand here — one node upstream of the Iceman's placement.
- **`G-Z6499` is proposed below `L166`**, on the branch toward the modern testers the block was
  built from. The Iceman is ancestral at it 10/0, but below-versus-parallel rests on parsimony alone.

Samples carrying all markers in an equivalence block cannot resolve their internal order. An
informative intermediate can split the block, including an ancient individual absent from a
modern-only comparison. This does not establish that such lineages are extinct, that no modern
intermediate exists, or that a commercial tree cannot incorporate ancient evidence.

## The ladder

```
G-PF3239                     I14677, I14678, I15942  (Sardinia)  ── stop here
  └── G-Z6219   provisional  MX210, SX10 (Oberbipp); Aesch12, Aesch13,
      │                      Aesch23 (Aesch); SNPRA61, SNPRA62 (Muttenz)
      └── G-L166             ÖTZI, I5118 (Hungary)
          ├── G-Z6208  refuted      ← the published label's node
          ├── G-Z6494              Iceman ancestral 3/3 — sibling branch
          │   └── G-Z6211          Iceman ancestral 2/2
          └── G-Z6499  provisional  Iceman ancestral 10/0 — modern branch
```

The topology is stated once, machine-readably, in
[`markers/tree_local.tsv`](../markers/tree_local.tsv), with a `status` column
(`published`/`provisional`/`putative`/`refuted`) that the tooling enforces: a terminal call may rest
only on a node's `defining` markers, and a `refuted` node can never appear on a path.

## Who the closest paternal relatives are

**`I5118` — Mezőcsát-Hörcsögös, Hungary, 3300–3000 BCE.** Near-contemporary with the Iceman.
Derived at two distinct transversion loci in the `L166` block (`L166` and `FGC5696`),
**one read each**, no conflicting read —
a `DERIVED_1read_transversion` call under rules fixed before the data existed. Backbone `P15`,
`PF3147`, `L91`, `PF3239`, `Z6043` all derived with no ancestral read. The call is a single molecule
per site and is reported as exactly that. In the local scorer, `L166` is defining and `FGC5696`
is an equivalent block marker. The tool also returns a thinner `L166` placement for `Aesch7`,
with a weak ancestral observation; its caveat is described below. Thus `I5118` is not literally the
only other tool-returned `L166` placement.

**One node upstream, at `G-Z6219`: seven men across three Swiss Neolithic sites** — Oberbipp
(`MX210`, `SX10`), Aesch (`Aesch12`, `Aesch13`, `Aesch23`) and Muttenz (`SNPRA61`, `SNPRA62`). Each
is derived at `Z6219` and, where callable, ancestral at the `L166` transversions. They are the
upstream comparison set in this analysis, not demonstrated close genealogical relatives or
independent replications of the branch.

These are **not demonstrably independent patrilines** and this design cannot make them so. Kin
groups negative at third degree are not thereby Y-independent, and neither the source publication
nor the aggregator measures patrilineal sharing beyond that. The honest count of independent
lineages carrying the `Z6219`-derived / `L166`-ancestral state is **1**, observed at three sites.

## Who is ruled out, and why

| excluded | evidence |
|---|---|
| **Three Sardinians** — `I14677`, `I14678`, `I15942` | `I14677` and `I14678` place at `G-PF3239`, `I15942` at `G-L91` on coverage; all three are excluded from both `G-Z6219` and `G-L166`. 42 ancestral reads across three tested `L166`-block transversion loci, zero derived. `I14677` alone carries 11 reads at `L166`. |
| **Upstream Oberbipp/Aesch/Muttenz comparison set** | `Z6219`-derived / `L166`-ancestral observations recur across the three sites. Cohort totals are 59 ancestral reads at `L166` and zero derived, and 25 derived at `Z6219` with zero ancestral. Not every man is callable at both sites; see the individual placements and the Aesch7 caveat. |
| **Modern `G-L166` block comparison** | The Iceman is ancestral at `Z6499` 10/0, contrasting with the derived block used for comparison. The local tree's below-versus-parallel interpretation remains provisional. |
| **`CGG017683`** (Crimea) | No placement power in the public deposit examined — zero reads at `Z6219`, `L166`, `L167` and `Z6499`, verified with all filters off. Of ~20 Iceman novel sites it covers 3, ancestral at all 3. Revisit only for a genuinely new input or correction. |
| **`UNTA58_68Sk1`** (Augsburg) | No power. Zero reads across the nine-marker `L166` panel tested. |

The usual C>T/G>A deamination mechanism cannot explain these transversion states at `L166` (`C>A`)
or `L167` (`T>A`), nor manufacture the ancestral C at `Z6499` (`C>T`). This is a directional damage
check, not a blanket exclusion of sequencing error, contamination or mismapping.

## What is provisional, and what would overturn it

- **`G-Z6219`** is supported at three sites and 25 derived reads with zero ancestral, and was
  pre-registered before the Aesch/Muttenz data was read. It rests on **one** demonstrably
  independent lineage. A single individual `Z6219`-ancestral and `L166`-derived would refute it.
- **`G-Z6499`'s direction** — below `L166` rather than parallel to it — **rests on parsimony alone.**
  No sample carrying `Z6499`-derived and `L166`-ancestral has been looked for. Its relation to
  `G-Z6494` is unresolved.
- **`Z6135` and `Z6209`** may belong beside `Z6219`. That is a post-hoc scan result over 22
  positions, unregistered, and carries no weight; they are recorded as `equivalent` markers, which
  the scorer counts but never lets deepen a placement.
- **`Aesch7`** places at `G-L166` off one derived `L167` read while carrying one low-power ancestral
  read at `L166`. Both call rules are correct and together they produce a placement from a single
  molecule. It prints `weak_ancestral=1` and is the only such sample in either cohort.

## Where the published labels went wrong

Three occurrences, three different columns, one mechanism: **the qualifier and the evidence live in
a different field from the label, and only the label propagates.**

The `haplogroup.info` "All Ancient DNA" compilation carries, per row, a `SNP-positive` list of every
derived call. For these samples that field is *correct* and disagrees with the row's own label:

| sample | deepest node in `SNP-positive` | label field | read-level verdict |
|---|---|---|---|
| `I14677`, `I14678`, `I15942` | `G2a2a1a2a1` (`PF3239`) | `G2a2a1a2a1a` | ancestral at `L166` |
| `I5118` | `G2a2a1a2a1a` incl. `L166:23989884C->A` | `G2a2a1a2a1a` | **derived** at `L166` |
| `E09538` | `FGC5671`, `Z6219` — *not* `L166` | `G2a2a1a2a1a` | the intermediate state |

The three samples the reads falsify are exactly the three whose own evidence field stops at
`PF3239`. **Rule: read `SNP-positive`, ignore `Y-Haplotree-Variant`.**

This is not attributable to any publication. `Source` names the paper the *sample* came from, not
the origin of the Y call; the compilation states its own method. Where the chain *has* been traced,
the over-calling has been downstream of the publications every time — Furtwängler 2020 assigns the
Oberbipp males `PF3239` and the aggregator moved them deeper. The Iceman's own `G-Z6208` shows
publications are not immune.

## Reproducing this

Everything above regenerates from committed tables. Run from the repository root:

```bash
python3 annotate/test_ytree.py                      # 40 checks on the tree file and scorer

python3 annotate/y_tree_place.py --sample Iceman \
  --pileup iceman-y/results/iceman_y_L166_evidence.tsv \
           iceman-y/results/z6219_node/iceman_yfull_L166_defining.tsv \
  --out iceman-y/results/placement/iceman_placement.tsv

python3 annotate/y_tree_place.py --pileup iceman-y/results/swiss15/swiss_*.tsv \
  --out iceman-y/results/placement/oberbipp15_placement.tsv

python3 annotate/y_tree_place.py --pileup iceman-y/results/testC/testC_*.tsv \
  --out iceman-y/results/placement/aesch_muttenz_placement.tsv

python3 annotate/y_tree_place.py \
  --pileup iceman-y/results/unhedged/unhedged_*.tsv \
           iceman-y/results/testB_unhedged/unhedged_*.tsv \
  --out iceman-y/results/placement/unhedged_placement.tsv

python3 annotate/ytree.py --newick iceman-y/results/placement/tree_local.nwk \
  --markers-out iceman-y/results/placement/tree_local_markers.tsv
```

`*_placement.tsv.summary.tsv` is one row per sample: terminal node, its status, what the sample is
excluded from, any caveats, and the full path. `tree_local.nwk` is the tree in Newick for external
placement tools (pathPhynder/phynder), refuted nodes dropped.

## The paperwork

The pre-registrations are the record that each prediction was written down before its data was
read; commit timestamps are what makes that claim checkable. They are not required reading.

- [`prereg/Z6219_node.md`](prereg/Z6219_node.md) — the `Z6219` node, plus extensions E1–E3
  (`FGC5671`/`Z6499` in the Iceman; the Crimean sample; the uniqueness filter)
- [`prereg/testC_aesch_muttenz.md`](prereg/testC_aesch_muttenz.md) — Test C, the Aesch/Muttenz
  cohort that took `Z6219` from one site to three
- [`prereg/swiss_neolithic_L166.md`](prereg/swiss_neolithic_L166.md) — the Oberbipp series
- [`prereg/unhedged_L166.md`](prereg/unhedged_L166.md) — the Sardinians and `I5118`
- [`prereg/uniqueness_filter.md`](prereg/uniqueness_filter.md) — a filter tested and rejected
- [`PROTOCOL_extending_analyses.md`](PROTOCOL_extending_analyses.md) — how a new question gets added

The dated [candidate assessment](research/otzi_literature_candidates.md),
[accession crosswalk](research/otzi_literature_candidates.tsv), and
[research priorities/source-watch note](research/priorities_and_source_watch.md) preserve the Work
assessments. They nominate prerequisites and possible follow-ups, not new placements. Their source
commits and the unverified watch-scheduling claim are recorded in [SITUATION.md](../SITUATION.md).

Working notes, gotchas and the full command history stay in the repository root
([`NOTES.md`](../NOTES.md), [`RUNLOG.md`](../RUNLOG.md), [`METHODS.md`](../METHODS.md)). Those files
predate this investigation and are shared with the mapping work; there is no clean line to split
them on, so they were left where they are.

## What this does not claim

Nothing here speaks to where anyone came from. The geographic spread — Tyrol, Hungary, three Swiss
sites, Sardinia — invites a migration reading and the pre-registrations refuse it in advance. These
are Y-chromosome placements on a tree, and a shared patriline is not a shared population.
