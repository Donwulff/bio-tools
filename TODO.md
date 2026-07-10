# TODO

1. `mapping/revert-bam.sh`: guard fastp usage. If fastp is enabled, require `--interleaved_out` or switch to paired FASTQ output. Otherwise skip fastp to preserve pairing.
2. `mapping/revert-bam.sh`: optional `NO_MARKDUP=1` path to skip MarkDuplicates when singleton-heavy data makes duplicate metrics unreliable.
3. `util/build_env.sh`: pass build flags explicitly on `make install` (use `sudo env CFLAGS=... LDFLAGS=... LIBS=...`). Avoid `sudo -E` which drops env and triggers rebuilds without zlib.
4. `util/build_env.sh`: standardize zlib linking across tools (bedtools/bwa/samtools/bcftools) when using `../zlib/libz.a`.
5. `util/build_env.sh`: Node build should pin a supported Python (`PYTHON=/usr/bin/python3.11` or `--python=...`) to avoid Python 3.13 configure failures.
