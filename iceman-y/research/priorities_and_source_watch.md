# Research value, candidate evidence and source watch

Decision note dated 2026-10-28. Baseline: `43273c02ce4a0bd788118d4da6c2a46460894d29`.

This records how to choose further work. It inherits the existing findings and their limits; it does not rerun analyses, change sample calls or tree placement, amend a preregistration, or close a registered test. The current literature/candidate search and discriminator search are separate workstreams; their outcomes have not been incorporated here.

## Decide value before checking whether a run is possible

Every proposed experiment must answer three separate questions:

1. **Value for the Iceman investigation:** what particular conclusion could change, under which observable outcomes?
2. **Value for reusable ancient-DNA tooling:** what general capability, failure mode or known-answer validation would it add?
3. **Feasibility:** which exact reads, source definitions, reference sequence and tool behavior are needed, and are they available?

An experiment may be feasible but not worth repeating. It may be useful for other ancient DNA while adding nothing to the Iceman case. Neither a new reference nor more computation recovers molecules already discarded upstream. Keep scientific value and input availability in separate fields rather than using either as a proxy for the other.

For each proposal record: question; value for Iceman; value for tooling; decisive observation; opposite outcome; no-power outcome; existing evidence; required inputs; expected effort; and decision (`run`, `prerequisite first`, `defer`, `not worth running for this aim`). A negative value assessment is not an automatic abandonment of a registered test: any abandonment must be explicit in its governing record.

## Current phase assessment

These are bounded planning judgments, not new experimental findings.

| Phase or purpose | Value for Iceman | Value for reusable tooling | Next decision |
|---|---|---|---|
| Remap the same selected reads solely to recover reads discarded by the original pipeline | None for that recovery objective: the rejected molecules are absent | A different, explicit reference-placement benchmark could still be useful | Do not rerun for the impossible recovery objective; name any genuinely different endpoint first |
| Repeat the established reference-recovery comparison under unchanged inputs/settings | Low; the comparison is already recorded | Reproduction is useful only for a concrete regression, implementation change or unresolved discrepancy | Reuse the committed evidence; do not relaunch a generic reference comparison |
| Off-target mapping onto chrY: native non-Y reads accepted at Y markers | Indirect and conditional: worthwhile if the result can change an interpretation of retained evidence or a proposed Y-only workflow | A useful candidate test of directional specificity beyond Y-origin recovery | Assess the decision it would change, then implement/register the bounded experiment if justified |
| Reconcile the authoritative Test B roster | Clarifies what was covered by the claim; new genotype work matters only if newly resolved positions can be informative | Improves source membership and missingness accounting | Resolve provenance first; do not assume a full rerun or new evidence follows |
| Inspect existing tables for the Z6499-derived / L166-ancestral discriminator | Potentially useful at low cost; a separate search is already in progress | No generic-tool dependency | Consume that search's result once; do not launch a duplicate |
| Find independent comparison samples | Potentially high if they provide a discriminating state with adequate coverage and defensible paternal independence | Can also supply independently documented validation cases | Rank actual individuals and available observations before staging reads |
| Path and history cleanup | Useful when it unblocks a named reproduction or handoff | Ordinary maintainability | Defer broad cleanup while research decisions remain open |

The [four-reference mappability record](../../NOTES.md#sensitivity-test-what-fraction-of-reads-can-even-come-back-2026-07-26) is a synthetic recovery comparison, not automatically a comparison of complete original-sample remappings. The [later correction](../../NOTES.md#the-chry-only-finding-measured-recovery-not-contamination-2026-07-27) limits what the chrY-only result establishes. The directional specificity experiment has not been run in this note.

The [read-selection assessment](../../NOTES.md#every-sample-in-this-project-came-from-a-mapped-only-deposit-2026-07-26) and [interpretation boundaries](../../METHODS.md#interpretation-boundaries) remain operative. The unresolved ten source groups in Test B, lack of reads at a resolved locus, reference incompatibility and QC rejection are different limitations. Do not infer one from another.

## Third-party assignments are candidate leads

The project has already investigated several claims associating ancient samples with the Iceman's Y branch. Future research should inherit the results below rather than make the user establish them again.

| Case | Inherited project finding | Rule for future candidate searches |
|---|---|---|
| `I14677`, `I14678`, `I15942` | The displayed L166 assignments fail the recorded read-level check: the three Sardinians provide 42 ancestral reads across the three tested L166-defining transversions, with no derived reads. The compilation's listed positive SNP evidence stops above L166. | Do not resurface the same assignments as positive L166 matches without concrete new evidence or a different dataset |
| `I5118` | The L166 state verifies under the standing single-read transversion rule: one derived read at L166 and one at FGC5696, without conflicting reads | Preserve the one-molecule-per-site limitation; it is a verified project observation, not a newly discovered candidate |
| Oberbipp/Aesch/Muttenz | Recorded Z6219-derived / L166-ancestral observations; additional sites do not demonstrate independent paternal lineages | Do not turn site count into independent replication or mark P3 met from it |
| `CGG017683` | No power at the decisive positions in the public deposit already examined | Revisit only for a named new input or evidence that changes that limit |
| `UNTA58_68Sk1` / `E09538` | One individual under two aggregation IDs; the staged UNTA input has no coverage at all nine L166-defining markers | Deduplicate identity; keep no-power distinct from a falsified genotype assignment |

Evidence is already in the [project verdict](../README.md#who-is-ruled-out-and-why), [unhedged-set assessment](../../NOTES.md#the-unhedged-g-l166-set-the-label-fails-in-both-directions-2026-07-26), [source-field comparison](../../NOTES.md#the-unhedged-g-l166-set-the-label-fails-in-both-directions-2026-07-26), and committed [unhedged results](../results/unhedged/). The Sardinian pair's registered observation count also remains unchanged.

The documented failure is specific: uncertainty or SNP evidence can be retained in one field while a stronger terminal label propagates from another. This does not establish intent, and it does not make every enthusiast result wrong.

- Preserve exact labels, provisional flags, aliases, `SNP-positive` evidence and source dates.
- A compilation's `Source` field may identify the sample publication rather than the source of its Y assignment. Do not attribute the curator's call to the publication without tracing it.
- Named SNP evidence is more informative than a displayed long-form label. An upstream positive call is not a terminal placement; a no-call is not an ancestral call.
- Use compilers and social posts to discover candidates, then trace an individual to primary evidence and actual available files.
- Revisit a ruled-out or no-power sample only when an explicit new input, correction or relevant new observation justifies it. A renamed row or repeated social claim is not new evidence.

## Comparison samples can serve different purposes

Keep candidate research in one workstream, but give every individual a declared role:

1. **Branch discrimination:** could its observed marker states distinguish the proposed placements?
2. **Pipeline verification:** does independently documented evidence give a known-answer check of mapping, filtering or placement? Discovery labels alone are not ground truth.
3. **Historical comparison:** is it interesting context without resolving the current hypotheses?

A sample useful for the second or third role need not be a close paternal relative. Keep those roles distinct when ranking cost against value. Broad G membership, culture or geography alone does not establish closeness.

## Periodic source watch

A monthly watch is configured for the first day of the month, using the user's Europe/Helsinki timezone. It checks:

- Official ISOGG haplogroup G material for additions, removals, aliases, changed parentage and provisional/reliability flags below or splitting the L91 / PF3239 / L166 path.
- The project's existing candidate-discovery source, [All Ancient DNA](https://haplogroup.info/all-ancient-dna.txt), and its [haplotree.info front-end](https://haplotree.info/maps/ancient_dna/), for genuinely new individuals or materially revised evidence/data releases.

Report exact source versions/dates and SNP/sample identities. A newly named downstream tree node is not evidence that the Iceman carries its defining allele; direct supporting evidence must be distinguished from nomenclature changes. Do not limit discovery to an exact L166 label when shallower labels may hide usable coverage.

Compare with existing project findings and earlier watch results. Establish a dated baseline on the first successful check. Notify only about substantive changes, or a persistent access blocker once until it changes; do not report access failure as evidence of no change. Record why each update matters for Iceman, reusable tooling, or neither yet.

This watch authorizes source checking and notification, not new genotyping, sequence downloads, repository edits or changes to recorded outcomes. The suggested literature/candidate research remains a separate bounded task. Existing [extension rules](../PROTOCOL_extending_analyses.md) and [preregistrations](../prereg/) govern any subsequent claim-bearing analysis.
