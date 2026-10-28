# TODO

Status reviewed 2026-10-28. See [SITUATION.md](SITUATION.md) for research/PR status and unresolved
registered work. Items below concern reusable tooling; they are not authorization to change
registered thresholds or rerun sample analyses.

1. Completed: `mapping/revert-bam.sh` includes fastp `--interleaved_out`. Keep the pairing regression
   checklist; this fix does not prove that mixed SE/PE chunk boundaries are safe.
2. Completed under a different option name: `ENABLE_MARKDUP=false` skips both duplicate markers;
   `ENABLE_ALIGN` also gates alignment. No additional `NO_MARKDUP` variable is needed.
3. `util/build_env.sh`: pass build flags explicitly on `make install` (use `sudo env CFLAGS=... LDFLAGS=... LIBS=...`). Avoid `sudo -E` which drops env and triggers rebuilds without zlib.
4. `util/build_env.sh`: standardize zlib linking across tools (bedtools/bwa/samtools/bcftools) when using `../zlib/libz.a`.
5. `util/build_env.sh`: Node build should pin a supported Python (`PYTHON=/usr/bin/python3.11` or `--python=...`) to avoid Python 3.13 configure failures.
6. Liftover: finish consolidating reference discovery, configuration, Picard/GATK execution and
   validation. Shared `lib_liftover.sh` exists, but the drivers still duplicate operations. Validate
   nonempty sidecars, matching dictionaries, staged compression/indexes and per-input output names.
7. gVCF: separate ordinary variant conversion from reference-block transformation. Picard block
   preservation is unvalidated; require a block-aware implementation and meaningful fixtures before
   advertising that capability. See [annotate/LIFTOVER.md](annotate/LIFTOVER.md).
8. `compare_marker_status_sets.py`: reject or explicitly reconcile duplicate identities/loci/alleles;
   add callable-on-both and expected-path denominators, allele depths, and lost-liftover accounting.
   Validate cross-assembly comparisons by verified identity, not coordinate equality or ID alone.
9. Marker resolution: interior gVCF lookup and all-state quality filtering are implemented with
   synthetic regressions (`annotate/test_y_haplo_from_markers.py`). Existing sample outputs were
   not regenerated. Remaining: exact panel hashes, broader reference/allele-orientation validation,
   and matched-panel expected-path QC; block-indel interpretation remains unsupported.
10. Ranking/reports: reduce stale partial reruns and output-schema drift; separate descriptive
    label/marker scores from topology-based placement. Do not present DP sums as supporting-molecule
    counts or tune scores on a known sample and report that as independent validation.
