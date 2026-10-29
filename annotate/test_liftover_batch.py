#!/usr/bin/env python3
"""Synthetic batch lifecycle tests; no private datasets or downloads."""

import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
import zlib

import liftover_to_hg38_batch as batch

try:
    import pysam
except ImportError:
    pysam = None


HEADER = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chrY,length=16>\n'
    '##contig=<ID=decoy,length=16>\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    '##ALT=<ID=NON_REF,Description="Unused declaration is allowed">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n'
)
ROW = 'chrY\t2\t.\tC\tG\t50\tPASS\t.\tGT\t1\n'

FAKE_JAVA = r'''
import gzip
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import pysam

opts = dict(arg.split("=", 1) for arg in sys.argv if "=" in arg)
source = Path(opts["I"])
out = Path(opts["O"])
reject = Path(opts["REJECT"])
with open(os.environ["TEST_CALLS"], "a") as calls:
    calls.write(source.name + "\n")
mode = os.environ.get("TEST_MODE", "success")
if mode == "sleep":
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    Path(os.environ["TEST_READY"]).write_text(f"{os.getpid()} {child.pid}")
    time.sleep(60)
with source.open("rb") as handle:
    zipped = handle.read(2) == b"\x1f\x8b"
opener = gzip.open if zipped else open
with opener(source, "rt") as handle:
    lines = handle.readlines()
header = "".join(line for line in lines if line.startswith("#"))
rows = [line for line in lines if not line.startswith("#")]
if mode == "wrong_samples":
    header = header.replace("SYNTHETIC", "CHANGED")
if mode == "lost":
    rows = rows[:-1]
if mode == "doubled":
    rows = rows * 2
if mode == "block":
    rows = [row.replace("\tG\t", "\t<NON_REF>\t") for row in rows]
with pysam.BGZFile(str(out), "w") as handle:
    handle.write((header + ("" if mode == "all_rejected" else "".join(rows))).encode())
with pysam.BGZFile(str(reject), "w") as handle:
    handle.write((header + ("".join(rows) if mode == "all_rejected" else "")).encode())
if mode == "mutate_input":
    with source.open("a") as handle:
        handle.write("\n")
if mode == "populate_final":
    target = Path(os.environ["TEST_FINAL"])
    target.mkdir(exist_ok=True)
    (target / "keep").write_text("not yours")
if mode == "corrupt":
    out.write_bytes(out.read_bytes()[:-12])
if mode == "fail" or (mode == "fail_second" and source.name.startswith("second")):
    print("synthetic backend failure", flush=True)
    sys.exit(23)
'''


class VcfScanTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="liftover-scan-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_content_guard_independent_of_name_and_compression(self):
        for alt, info in (("<NON_REF>", "."), ("<*>", "END=3"),
                          ("G,<NON_REF>", "."), (".", "END=4"), (".", "END=2")):
            for compressed in (False, True):
                with self.subTest(alt=alt, info=info, compressed=compressed):
                    path = self.root / "ordinary.vcf"
                    text = HEADER + ROW.replace("\tG\t", f"\t{alt}\t").replace(
                        "\tPASS\t.\t", f"\tPASS\t{info}\t")
                    path.write_bytes(gzip.compress(text.encode()) if compressed else text.encode())
                    with self.assertRaisesRegex(batch.BatchError, "unsupported gVCF"):
                        batch.scan_vcf(path)

    def test_unused_gvcf_header_spanning_deletion_and_symbolic_variant_allowed(self):
        for alt, info in (("G", "."), ("*", "."), ("<DEL>", "END=4")):
            with self.subTest(alt=alt):
                path = self.root / "calls.gvcf.gz"
                path.write_text(HEADER + ROW.replace("\tG\t", f"\t{alt}\t").replace(
                    "\tPASS\t.\t", f"\tPASS\t{info}\t"))
                self.assertEqual(batch.scan_vcf(path)["records"], 1)

    def test_truncated_gzip_rejected(self):
        path = self.root / "ordinary.vcf"
        good = gzip.compress((HEADER + ROW).encode())
        for data in (good[:-8], good[:-8] + b"\0" * 8):
            with self.subTest(trailer=data[-8:]):
                path.write_bytes(data)
                with self.assertRaises((EOFError, OSError, zlib.error)):
                    batch.scan_vcf(path)

    def test_empty_and_malformed_headers_records_rejected(self):
        for text in ("", "not a vcf\n", HEADER + "chrY\t2\n", HEADER + HEADER,
                     HEADER.replace("\tSYNTHETIC", "\tSYNTHETIC\tSYNTHETIC"),
                     HEADER + ROW.replace("\t2\t", "\t0\t"),
                     HEADER + ROW.replace("\tPASS\t.\t", "\tPASS\tEND=1\t")):
            with self.subTest(text=text[-50:]):
                path = self.root / "bad.vcf"
                path.write_text(text)
                with self.assertRaises(batch.BatchError):
                    batch.scan_vcf(path)

    def test_site_only_and_header_only_inputs(self):
        path = self.root / "sites.vcf"
        path.write_text(HEADER.rsplit("\tFORMAT", 1)[0] + "\n")
        self.assertEqual(batch.scan_vcf(path), {"records": 0, "samples": []})

    def test_stems_strip_one_suffix_preserving_dots(self):
        for name in ("sample.part.vcf.gz", "sample.part.gvcf", "sample.part.g.vcf.gz"):
            self.assertEqual(batch.output_stem(Path(name)), "sample.part")
        with self.assertRaises(batch.BatchError):
            batch.output_stem(Path(".vcf"))


@unittest.skipUnless(pysam and all(shutil.which(t) for t in ("bcftools", "tabix", "samtools")),
                     "requires pysam, bcftools, tabix and samtools")
class BatchLifecycleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="liftover-batch-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        fake = self.bin / "java"
        fake.write_text(f"#!{sys.executable}\n" + FAKE_JAVA)
        fake.chmod(0o755)
        self.calls = self.root / "calls"
        self.ready = self.root / "ready"
        self.out = self.root / "final output"
        self.ref = self.root / "reference.fa"
        self.ref.write_text(">chrY\nACGTACGTACGTACGT\n>decoy\nACGTACGTACGTACGT\n")
        pysam.faidx(str(self.ref))
        self.ref.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chrY\tLN:16\n@SQ\tSN:decoy\tLN:16\n")
        self.jar = self.root / "picard.jar"
        self.jar.write_text("synthetic test backend")
        self.chain = self.root / "identity.chain"
        self.chain.write_text("synthetic test chain")
        self.sample = self.root / "sample.part with spaces.vcf"
        self.sample.write_text(HEADER + ROW)
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ["PATH"],
                        TEST_CALLS=str(self.calls), TEST_READY=str(self.ready), TEST_FINAL=str(self.out))
        # Ignore user-level resource/JVM settings in synthetic runs.
        self.env.pop("PICARD_JAVA_OPTS", None)

    def command(self, samples=None, extra=()):
        command = ["sh", str(batch.SCRIPT_DIR / "liftover_to_hg38_batch.sh"),
                   "--chain-hs1-to-hg38", str(self.chain), "--ref-hg38", str(self.ref),
                   "--picard-jar", str(self.jar), "--out-dir", str(self.out), *extra]
        for sample in samples if samples is not None else [self.sample]:
            command += ["--sample", str(sample)]
        return command

    def run_batch(self, samples=None, mode="success", extra=()):
        return subprocess.run(self.command(samples, extra), env=dict(self.env, TEST_MODE=mode),
                              text=True, capture_output=True, timeout=25)

    def stages(self):
        return list(self.root.glob(".final output.staging-*"))

    def assert_preflight_failed(self, result):
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertFalse(self.calls.exists())
        self.assertFalse(self.stages())
        self.assertFalse((self.root / ".final output.liftover.lock").exists())

    def test_success_publishes_complete_batch_and_preserves_dots_spaces(self):
        second = self.root / "second.gvcf.gz"
        second.write_bytes(gzip.compress((HEADER + ROW).encode()))
        result = self.run_batch([self.sample, second])
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads((self.out / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual(len(manifest["inputs"]), 2)
        self.assertFalse(self.stages())
        self.assertFalse((self.root / ".final output.liftover.lock").exists())
        with (self.out / "liftover_report.tsv").open() as handle:
            report = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(report), 2)
        for row in report:
            self.assertEqual(row["input_type"], "vcf")
            self.assertEqual(row["lifted_records"], "1")
            self.assertEqual(row["rejected_records"], "0")
            self.assertEqual(row["primary_records"], "NA")
            self.assertTrue(Path(row["lifted_vcf"]).is_file())
            self.assertTrue(Path(row["lifted_vcf"] + ".tbi").is_file())

    def test_existing_empty_directory_supported(self):
        self.out.mkdir()
        result = self.run_batch()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads((self.out / "manifest.json").read_text())["status"], "complete")

    def test_primary_subset(self):
        self.sample.write_text(HEADER + ROW + ROW.replace("chrY", "decoy"))
        result = self.run_batch(extra=["--primary-only"])
        self.assertEqual(result.returncode, 0, result.stderr)
        plan = json.loads((self.out / "manifest.json").read_text())["inputs"][0]
        self.assertEqual((plan["lifted_records"], plan["primary_records"]), (2, 1))

    def test_all_rejected_and_empty_vcf_are_valid(self):
        for empty in (False, True):
            with self.subTest(empty=empty):
                self.out = self.root / f"result-{empty}"
                self.sample.write_text(HEADER if empty else HEADER + ROW)
                result = self.run_batch(mode="all_rejected", extra=["--primary-only"])
                self.assertEqual(result.returncode, 0, result.stderr)
                plan = json.loads((self.out / "manifest.json").read_text())["inputs"][0]
                self.assertEqual(plan["lifted_records"], 0)
                self.assertEqual(plan["rejected_records"], 0 if empty else 1)
                self.assertEqual(plan["primary_records"], 0)

    def test_missing_later_input_does_not_touch_outputs_or_reference_metadata(self):
        Path(str(self.ref) + ".fai").unlink()
        result = self.run_batch([self.sample, self.root / "missing.vcf"])
        self.assert_preflight_failed(result)
        self.assertFalse(Path(str(self.ref) + ".fai").exists())
        self.assertFalse(self.out.exists())

    def test_malformed_later_input(self):
        second = self.root / "second.vcf"
        second.write_text(HEADER + ROW + "bad\n")
        self.assert_preflight_failed(self.run_batch([self.sample, second]))

    def test_reference_blocks_rejected_before_reference_writes(self):
        Path(str(self.ref) + ".fai").unlink()
        self.sample.write_bytes(gzip.compress((HEADER + ROW.replace("\tG\t", "\t<*>\t")).encode()))
        result = self.run_batch()
        self.assert_preflight_failed(result)
        self.assertIn("unsupported gVCF", result.stderr)
        self.assertFalse(Path(str(self.ref) + ".fai").exists())

    def test_input_basename_collision(self):
        other = self.root / "another"
        other.mkdir()
        second = other / self.sample.name
        second.write_text(HEADER + ROW)
        result = self.run_batch([self.sample, second])
        self.assert_preflight_failed(result)
        self.assertIn("colliding output basename", result.stderr)

    def test_input_aliases_and_hardlinks(self):
        for hardlink in (False, True):
            with self.subTest(hardlink=hardlink):
                alias = self.root / f"alias-{hardlink}.vcf"
                if hardlink:
                    os.link(self.sample, alias)
                else:
                    alias.symlink_to(self.sample)
                result = self.run_batch([self.sample, alias])
                self.assert_preflight_failed(result)
                self.assertIn("duplicate input", result.stderr)

    def test_populated_directory_preserves_every_file(self):
        self.out.mkdir()
        for name in ("liftover_report.tsv", "other.txt", "orphan.hg38.primary.vcf.gz.tbi"):
            (self.out / name).write_text("keep unchanged")
        self.assert_preflight_failed(self.run_batch(extra=["--primary-only"]))
        self.assertTrue(all(p.read_text() == "keep unchanged" for p in self.out.iterdir()))

    def test_symlink_destination_rejected(self):
        self.out.symlink_to(self.root / "absent", target_is_directory=True)
        self.assert_preflight_failed(self.run_batch())
        self.assertTrue(self.out.is_symlink())

    def test_existing_lock_preserved(self):
        lock = self.root / ".final output.liftover.lock"
        lock.write_text("existing owner")
        result = self.run_batch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("run lock exists", result.stderr)
        self.assertEqual(lock.read_text(), "existing owner")
        self.assertFalse(self.stages())

    def test_invalid_reference_metadata_retained_without_publication(self):
        dictionary = self.ref.with_suffix(".dict")
        dictionary.write_text("invalid existing metadata\n")
        result = self.run_batch()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())
        self.assertFalse(self.out.exists())
        self.assertEqual(dictionary.read_text(), "invalid existing metadata\n")
        stage = self.stages()[0]
        self.assertEqual(json.loads((stage / "manifest.json").read_text())["status"], "failed")
        self.assertTrue((stage / "logs/reference.log").is_file())

    def test_missing_reference_index_created_only_after_valid_preflight(self):
        Path(str(self.ref) + ".fai").unlink()
        result = self.run_batch()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(Path(str(self.ref) + ".fai").is_file())

    def test_invalid_java_options_fail_before_writes(self):
        result = self.run_batch(extra=["--java-opts", '"unterminated'])
        self.assert_preflight_failed(result)
        self.assertIn("No closing quotation", result.stderr)

    def test_index_failure_cannot_publish(self):
        tabix = self.bin / "tabix"
        tabix.write_text("#!/bin/sh\nexit 19\n")
        tabix.chmod(0o755)
        result = self.run_batch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("command exited 19", result.stderr)
        self.assertFalse(self.out.exists())
        self.assertEqual(json.loads((self.stages()[0] / "manifest.json").read_text())["status"], "failed")

    def test_backend_failures_never_publish_and_keep_logs_manifest(self):
        for mode, expected in (("fail", "command exited 23"), ("lost", "conservation"),
                               ("doubled", "conservation"), ("wrong_samples", "sample names"),
                               ("block", "unsupported gVCF"), ("corrupt", "ERROR:"),
                               ("mutate_input", "input changed")):
            with self.subTest(mode=mode):
                self.sample.write_text(HEADER + ROW)
                result = self.run_batch(mode=mode)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected, result.stderr)
                self.assertFalse(self.out.exists())
                stage = self.stages()[0]
                self.assertEqual(json.loads((stage / "manifest.json").read_text())["status"], "failed")
                self.assertTrue(list((stage / "logs").glob("*.tools.log")))
                self.assertFalse((self.root / ".final output.liftover.lock").exists())
                shutil.rmtree(stage)

    def test_second_backend_failure_keeps_first_unpublished(self):
        second = self.root / "second.vcf"
        second.write_text(HEADER + ROW)
        result = self.run_batch([self.sample, second], mode="fail_second")
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.out.exists())
        manifest = json.loads((self.stages()[0] / "manifest.json").read_text())
        self.assertEqual(manifest["inputs"][0]["status"], "validated")
        self.assertEqual(manifest["status"], "failed")

    def test_concurrent_destination_content_not_overwritten(self):
        result = self.run_batch(mode="populate_final")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.out / "keep").read_text(), "not yours")
        self.assertEqual([p.name for p in self.out.iterdir()], ["keep"])
        self.assertTrue(self.stages())

    def test_sigterm_stops_tool_group_and_retains_stage(self):
        for signum in (signal.SIGTERM, signal.SIGINT):
            with self.subTest(signal=signum):
                self.ready.unlink(missing_ok=True)
                proc = subprocess.Popen(self.command(), env=dict(self.env, TEST_MODE="sleep"),
                                        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                try:
                    deadline = time.monotonic() + 15
                    while not self.ready.exists() and time.monotonic() < deadline and proc.poll() is None:
                        time.sleep(0.05)
                    self.assertTrue(self.ready.exists(), "backend did not start")
                    pids = [int(s) for s in self.ready.read_text().split()]
                    proc.send_signal(signum)
                    stdout, stderr = proc.communicate(timeout=10)
                    self.assertEqual(proc.returncode, 128 + signum, stdout + stderr)
                    for pid in pids:
                        stat = Path(f"/proc/{pid}/stat")
                        self.assertTrue(not stat.exists() or stat.read_text().split()[2] == "Z")
                    self.assertFalse(self.out.exists())
                    stage = self.stages()[0]
                    self.assertEqual(json.loads((stage / "manifest.json").read_text())["status"],
                                     "interrupted")
                    self.assertFalse((self.root / ".final output.liftover.lock").exists())
                    shutil.rmtree(stage)
                finally:
                    if proc.poll() is None:
                        proc.kill()
                    proc.communicate()


if __name__ == "__main__":
    unittest.main()
