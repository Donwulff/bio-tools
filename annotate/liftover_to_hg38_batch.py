#!/usr/bin/env python3
"""Preflight ordinary VCF batches, validate staged liftover, then publish together."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import zlib

from liftover_resources import ResourceError, config_environment, find_reference
from liftover_backend import DEFAULT_JAVA_OPTS, MODES, select_backend


SCRIPT_DIR = Path(__file__).resolve().parent
PRIMARY_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
VCF_SUFFIXES = (".g.vcf.gz", ".gvcf.gz", ".vcf.gz", ".g.vcf", ".gvcf", ".vcf")


class BatchError(Exception):
    pass


class Interrupted(Exception):
    def __init__(self, signum):
        self.signum = signum
        super().__init__(f"signal {signum}")


def require_file(path):
    path = Path(path).expanduser().resolve(strict=True)
    if not path.is_file() or not os.access(path, os.R_OK) or not path.stat().st_size:
        raise BatchError(f"not a readable, nonempty file: {path}")
    return path


def snapshot(path):
    st = path.stat()
    return [st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns]


def assert_unchanged(path, expected):
    if snapshot(path) != expected:
        raise BatchError(f"input changed during this run: {path}")


def scan_vcf(path):
    """Structural scan; Picard validates remaining VCF semantics inside staging.

    Sniff compression rather than trusting suffixes. Reading gzip to EOF also
    checks its checksum/trailer; no truncation tolerance is enabled.
    """
    with path.open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    columns = None
    count = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if lineno == 1 and not line.startswith("##fileformat=VCFv4."):
                raise BatchError(f"{path}: missing VCFv4 fileformat header")
            if line.startswith("##") and columns is None:
                continue
            if line.startswith("#CHROM\t") and columns is None:
                columns = line.split("\t")
                fixed = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                if columns[:8] != fixed or (len(columns) != 8 and
                        (len(columns) < 10 or columns[8] != "FORMAT")):
                    raise BatchError(f"{path}:{lineno}: malformed VCF columns")
                samples = columns[9:]
                if len(set(samples)) != len(samples) or any(not s for s in samples):
                    raise BatchError(f"{path}:{lineno}: empty or duplicate sample names")
                continue
            if columns is None or line.startswith("#"):
                raise BatchError(f"{path}:{lineno}: misplaced/missing VCF header")
            fields = line.split("\t")
            if len(fields) != len(columns):
                raise BatchError(f"{path}:{lineno}: record/header column count differs")
            try:
                pos = int(fields[1])
                if pos < 1 or not all(fields[i] for i in (0, 3, 4)):
                    raise ValueError()
                ends = [v[4:] for v in fields[7].split(";") if v.startswith("END=")]
                if len(ends) > 1 or (ends and int(ends[0]) < pos):
                    raise ValueError()
            except ValueError as exc:
                raise BatchError(f"{path}:{lineno}: invalid POS/REF/ALT/END") from exc
            if {"<NON_REF>", "<*>"}.intersection(fields[4].split(",")) or (
                    fields[4] == "." and ends):
                raise BatchError(
                    f"{path}:{lineno}: unsupported gVCF/reference-confidence record; "
                    "Picard does not provide validated block conversion. Use an ordinary "
                    "variant VCF, or query the native gVCF with lifted markers instead."
                )
            count += 1
    if columns is None:
        raise BatchError(f"{path}: missing #CHROM header")
    return {"records": count, "samples": columns[9:]}


def output_stem(path):
    name = path.name
    for suffix in VCF_SUFFIXES:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    if not name or any(ord(c) < 32 or ord(c) == 127 for c in name):
        raise BatchError(f"invalid output basename: {path}")
    return name


def check_destination(path):
    if path.is_symlink() or (path.exists() and
                            (not path.is_dir() or any(path.iterdir()))):
        raise BatchError(f"output directory must be new or empty (no overwrite): {path}")


def run_command(argv, log, env_updates=None):
    """Stop the entire child group on cancellation, including shell-helper children."""
    with log.open("ab") as stream:
        for key, value in (env_updates or {}).items():
            stream.write((f"{key}={shlex.quote(value)}\n").encode())
        stream.write(("$ " + shlex.join([str(a) for a in argv]) + "\n").encode())
        stream.flush()
        proc = subprocess.Popen(argv, stdout=stream, stderr=subprocess.STDOUT,
                                start_new_session=True, env={**os.environ, **(env_updates or {})})
        try:
            code = proc.wait()
        except BaseException:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
                proc.wait(timeout=5)
            except ProcessLookupError:
                pass
            except subprocess.TimeoutExpired:
                pass
            # The leader may exit before a descendant; clean up that group too.
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
            raise
    if code:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        raise BatchError(f"command exited {code}; see {log}")


def validate_output(path, expected_samples):
    import pysam

    info = scan_vcf(path)
    if info["samples"] != expected_samples:
        raise BatchError(f"output sample names/order changed: {path}")
    # Decode all records and read through the generated index before publication.
    with pysam.VariantFile(str(path)) as handle:
        decoded = sum(1 for _ in handle)
    if decoded != info["records"]:
        raise BatchError(f"record decoding/count mismatch: {path}")
    with pysam.TabixFile(str(path), index=str(path) + ".tbi") as index:
        indexed = sum(1 for contig in index.contigs for _ in index.fetch(contig))
    if indexed != decoded:
        raise BatchError(f"index/data count mismatch: {path}")
    return decoded


def write_manifest(stage, manifest):
    tmp = stage / "manifest.json.tmp"
    with tmp.open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    tmp.replace(stage / "manifest.json")


def parser():
    p = argparse.ArgumentParser(description=__doc__, epilog=(
        "Full gVCF/reference-block conversion is unsupported, regardless of filename. "
        "Failed runs retain staging/logs; the final directory stays untouched."))
    p.add_argument("--chain-hs1-to-hg38", required=True)
    p.add_argument("--ref-hg38")
    p.add_argument("--ref-root", default=os.environ.get("REF_ROOT") or
                   os.environ.get("BIO_TOOLS_INDEX_DIR") or os.environ.get("INDEX_DIR"))
    p.add_argument("--allow-extended-ref", action="store_true",
                   help="explicitly permit extended/decoy-enriched reference filenames")
    p.add_argument("--picard-jar", default=os.environ.get("PICARD_JAR"))
    p.add_argument("--backend", choices=MODES, default=os.environ.get("LIFTOVER_BACKEND", "auto"))
    p.add_argument("--java-opts", default=os.environ.get("PICARD_JAVA_OPTS", DEFAULT_JAVA_OPTS))
    p.add_argument("--out-dir")
    p.add_argument("--primary-only", action="store_true")
    p.add_argument("--sample", action="append", required=True, help="ordinary VCF; repeat for a batch")
    return p


def preflight(args):
    for command in ("bcftools", "tabix", "samtools", "sh"):
        if not shutil.which(command):
            raise BatchError(f"missing command: {command}")
    try:
        import pysam  # noqa: F401 - dependency preflight before any writes
    except ImportError as exc:
        raise BatchError("pysam is required (python3 -m pip install pysam)") from exc
    chain = require_file(args.chain_hs1_to_hg38)
    ref = require_file(find_reference("hg38", args.ref_hg38, args.ref_root,
                                     allow_extended=args.allow_extended_ref))
    backend = select_backend(args.backend, args.picard_jar, args.java_opts)
    destination = Path(args.out_dir or Path.cwd() / "experiments" /
                       datetime.now().strftime("liftover_hg38_%Y%m%d_%H%M%S")).expanduser().absolute()
    # Resolve parents but not the leaf: a destination symlink must not be followed.
    destination = destination.parent.resolve() / destination.name
    check_destination(destination)
    plans = []
    identities = set()
    names = set()
    for value in args.sample:
        original = Path(value).expanduser()
        path = require_file(original)
        state = snapshot(path)
        if tuple(state[:2]) in identities:
            raise BatchError(f"duplicate input (including aliases/hardlinks): {path}")
        identities.add(tuple(state[:2]))
        stem = output_stem(original)
        if stem in names:
            raise BatchError(f"colliding output basename: {stem}")
        names.add(stem)
        info = scan_vcf(path)
        assert_unchanged(path, state)
        plans.append({"input": str(path), "snapshot": state, "stem": stem, **info,
                      "lifted_vcf": f"{stem}.hg38.vcf.gz",
                      "reject_vcf": f"{stem}.hg38.reject.vcf.gz",
                      "primary_vcf": f"{stem}.hg38.primary.vcf.gz" if args.primary_only else None})
    return chain, ref, backend, destination, plans


def execute(args):
    chain, ref, backend, destination, plans = preflight(args)
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.parent / f".{destination.name}.liftover.lock"
    stage = None
    manifest = {"status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
                "output_directory": str(destination), "reference": str(ref), "chain": str(chain),
                "picard_jar": backend.jar, "java_options": list(backend.java_options),
                "backend": backend.metadata(), "inputs": plans}
    try:
        lock_fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise BatchError(f"run lock exists; inspect before removing: {lock}") from exc
    lock_identity = os.fstat(lock_fd)
    try:
        with os.fdopen(lock_fd, "w") as handle:
            handle.write(f"pid={os.getpid()}\noutput={destination}\n")
        check_destination(destination)
        stage = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent))
        (stage / "logs").mkdir()
        write_manifest(stage, manifest)
        print(f"Staging batch in {stage}", flush=True)
        # Keep one reference validator/preparer shared with the other drivers.
        run_command(["sh", "-eu", "-c",
                     '. "$1/lib_liftover.sh"; SCRIPT_DIR="$1"; ensure_ref_sidecars "$2" "$3" "$4" "$5"',
                     "liftover-reference", str(SCRIPT_DIR), str(ref), backend.jar or "", args.java_opts,
                     json.dumps(backend.metadata())],
                    stage / "logs/reference.log")
        for plan in plans:
            source = Path(plan["input"])
            assert_unchanged(source, plan["snapshot"])
            lifted = stage / plan["lifted_vcf"]
            rejected = stage / plan["reject_vcf"]
            log = stage / "logs" / f"{plan['stem']}.tools.log"
            print(f"Liftover: {source.name}", flush=True)
            argv, updates = backend.invocation("LiftoverVcf", ["--INPUT", source, "--OUTPUT", lifted,
                "--CHAIN", chain, "--REFERENCE_SEQUENCE", ref, "--REJECT", rejected,
                "--RECOVER_SWAPPED_REF_ALT", "true", "--WARN_ON_MISSING_CONTIG", "true"])
            run_command(argv, log, updates)
            for path in (lifted, rejected):
                run_command(["tabix", "-p", "vcf", str(path)], log)
            n_lifted = validate_output(lifted, plan["samples"])
            n_rejected = validate_output(rejected, plan["samples"])
            if n_lifted + n_rejected != plan["records"]:
                raise BatchError(f"record conservation failed for {source.name}: "
                                 f"{n_lifted} + {n_rejected} != {plan['records']}")
            n_primary = None
            if plan["primary_vcf"]:
                primary = stage / plan["primary_vcf"]
                run_command(["bcftools", "view", "-r", ",".join(PRIMARY_CONTIGS),
                             "-Oz", "-o", str(primary), str(lifted)], log)
                run_command(["tabix", "-p", "vcf", str(primary)], log)
                n_primary = validate_output(primary, plan["samples"])
                if n_primary > n_lifted:
                    raise BatchError(f"primary subset count exceeds lifted count: {primary}")
            assert_unchanged(source, plan["snapshot"])
            plan.update(lifted_records=n_lifted, rejected_records=n_rejected,
                        primary_records=n_primary, status="validated")
            write_manifest(stage, manifest)
        with (stage / "liftover_report.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", "input_type", "lifted_vcf", "reject_vcf", "lifted_records",
                             "rejected_records", "primary_records"])
            for plan in plans:
                writer.writerow([plan["input"], "vcf", str(destination / plan["lifted_vcf"]),
                                 str(destination / plan["reject_vcf"]), plan["lifted_records"],
                                 plan["rejected_records"], plan["primary_records"]
                                 if plan["primary_records"] is not None else "NA"])
        for plan in plans:
            assert_unchanged(Path(plan["input"]), plan["snapshot"])
        check_destination(destination)
        manifest.update(status="complete", finished_utc=datetime.now(timezone.utc).isoformat())
        write_manifest(stage, manifest)
        # Same-parent rename publishes all products together; nonempty targets cannot be replaced.
        stage.replace(destination)
        stage = None
        print(f"Done. Report: {destination / 'liftover_report.tsv'}")
    except BaseException as exc:
        if stage is not None:
            manifest.update(status="interrupted" if isinstance(exc, (Interrupted, KeyboardInterrupt))
                            else "failed", error=str(exc),
                            finished_utc=datetime.now(timezone.utc).isoformat())
            try:
                write_manifest(stage, manifest)
            except OSError as write_error:
                print(f"Could not update failure manifest: {write_error}", file=sys.stderr)
            print(f"Unpublished staging retained: {stage}", file=sys.stderr)
        raise
    finally:
        # Never remove a lock replaced by another process or administrator.
        try:
            current = lock.lstat()
        except FileNotFoundError:
            current = None
        if current and (current.st_dev, current.st_ino) == (lock_identity.st_dev, lock_identity.st_ino):
            lock.unlink()


def main():
    try:
        environment = config_environment()
        os.environ.clear()
        os.environ.update(environment)
    except (ResourceError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    args = parser().parse_args()
    def interrupt(signum, _frame):
        raise Interrupted(signum)
    signal.signal(signal.SIGTERM, interrupt)
    try:
        execute(args)
    except Interrupted as exc:
        print(f"ERROR: interrupted by signal {exc.signum}", file=sys.stderr)
        return 128 + exc.signum
    except KeyboardInterrupt:
        print("ERROR: interrupted", file=sys.stderr)
        return 130
    except (BatchError, OSError, EOFError, ValueError, zlib.error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
