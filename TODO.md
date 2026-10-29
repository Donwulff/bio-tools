# TODO

Status updated 2026-10-29. See [SITUATION.md](SITUATION.md) for research/PR status and unresolved
registered work. Items below concern reusable tooling; they are not authorization to change
registered thresholds or rerun sample analyses.

1. Completed: `mapping/revert-bam.sh` includes fastp `--interleaved_out`. Keep the pairing regression
   checklist; this fix does not prove that mixed SE/PE chunk boundaries are safe.
2. Completed under a different option name: `ENABLE_MARKDUP=false` skips both duplicate markers;
   `ENABLE_ALIGN` also gates alignment. No additional `NO_MARKDUP` variable is needed.
3. `util/build_env.sh`: pass build flags explicitly on `make install` (use `sudo env CFLAGS=... LDFLAGS=... LIBS=...`). Avoid `sudo -E` which drops env and triggers rebuilds without zlib.
4. `util/build_env.sh`: standardize zlib linking across tools (bedtools/bwa/samtools/bcftools) when using `../zlib/libz.a`.
5. `util/build_env.sh`: Node build should pin a supported Python (`PYTHON=/usr/bin/python3.11` or `--python=...`) to avoid Python 3.13 configure failures.
6. Liftover: shared reference/index/dictionary validation, safe reference staging, reference discovery
   and local configuration loading are implemented with synthetic tests. Remaining consolidation:
   backend discovery/execution still differs among drivers; unrelated workflows do not yet share a
   repository-wide resource registry. Batch preflight, content-based confidence-record rejection,
   staged validation and atomic directory publication are implemented with synthetic lifecycle tests.
   These publication guarantees do not yet cover other experiment drivers.
7. gVCF: ordinary batch conversion now refuses confidence records based on content. Picard block
   preservation remains unvalidated; require a block-aware implementation and meaningful fixtures before
   advertising that capability. See [annotate/LIFTOVER.md](annotate/LIFTOVER.md).
8. Comparator: duplicate/conflict safeguards, callable-only metrics and explicit expected-panel
   denominators are implemented and synthetic-tested. Remaining: allele-depth export (not in the
   status schema), trusted coordinate/allele crosswalks and exact panel hashes, and actual
   liftover-reject provenance. Missing table IDs must not be labelled rejected loci automatically.
9. Marker resolution: interior gVCF lookup and all-state quality filtering are implemented with
   synthetic regressions (`annotate/test_y_haplo_from_markers.py`). Existing sample outputs were
   not regenerated. Remaining: exact panel hashes, broader reference/allele-orientation validation,
   and matched-panel expected-path QC; block-indel interpretation remains unsupported.
10. Ranking/reports: comparator output/stale-report guards are implemented; other drivers still
    need stronger rerun freshness and schema tracking. Separate descriptive
    label/marker scores from topology-based placement. Do not present DP sums as supporting-molecule
    counts or tune scores on a known sample and report that as independent validation.
