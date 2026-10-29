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
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import pysam

tool = next(t for t in ("LiftoverVcf", "CreateSequenceDictionary") if t in sys.argv)
args = sys.argv[sys.argv.index(tool) + 1:]
opts = {args[i].removeprefix("--"): args[i + 1] for i in range(0, len(args), 2)}
if os.environ.get("TEST_BACKEND_ARGS"):
    with open(os.environ["TEST_BACKEND_ARGS"], "a") as log:
        log.write(json.dumps({"program": Path(sys.argv[0]).name, "argv": sys.argv[1:],
                              "JAVA_TOOL_OPTIONS": os.environ.get("JAVA_TOOL_OPTIONS")}) + "\n")
if tool == "CreateSequenceDictionary":
    if os.environ.get("TEST_MODE") == "fail_dict":
        sys.exit(23)
    subprocess.run([os.environ["TEST_SAMTOOLS"], "dict", "-o", opts["OUTPUT"], opts["REFERENCE"]], check=True)
    sys.exit(0)
source = Path(opts["INPUT"])
out = Path(opts["OUTPUT"])
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
                        TEST_CALLS=str(self.calls), TEST_READY=str(self.ready), TEST_FINAL=str(self.out),
                        TEST_SAMTOOLS=shutil.which("samtools"))
        # Ignore user-level resource/JVM settings in synthetic runs.
        self.env.pop("PICARD_JAVA_OPTS", None)
        self.env.pop("BIO_TOOLS_CFG", None)
        self.env.pop("LIFTOVER_BACKEND", None)

    def install_fake_wrappers(self):
        for name in ("gatk", "picard"):
            shutil.copy2(self.bin / "java", self.bin / name)

    def command(self, samples=None, extra=()):
        command = ["sh", str(batch.SCRIPT_DIR / "liftover_to_hg38_batch.sh"),
                   "--chain-hs1-to-hg38", str(self.chain), "--ref-hg38", str(self.ref),
                   "--picard-jar", str(self.jar), "--out-dir", str(self.out), *extra]
        for sample in samples if samples is not None else [self.sample]:
            command += ["--sample", str(sample)]
        return command

    def run_batch(self, samples=None, mode="success", extra=()):
        return subprocess.run(self.command(samples, extra), env=dict(self.env, TEST_MODE=mode),
                              cwd=self.root, text=True, capture_output=True, timeout=25)

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

    def test_cli_reference_and_backend_override_local_config(self):
        (self.root / "bio-tools.cfg").write_text(
            "REF_HG38=missing.fa\nPICARD_JAR=missing.jar\nPICARD_JAVA_OPTS=invalid\nLIFTOVER_BACKEND=gatk\n")
        result = self.run_batch(extra=["--java-opts", "-Xms1g -Xmx2g", "--backend", "picard-jar"])
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads((self.out / "manifest.json").read_text())
        self.assertEqual(manifest["reference"], str(self.ref))
        self.assertEqual(manifest["java_options"], ["-Xms1g", "-Xmx2g"])

    def test_extended_reference_requires_explicit_opt_in(self):
        extended = self.root / "custom.H3630.fa"
        extended.symlink_to(self.ref)
        result = self.run_batch(extra=["--ref-hg38", str(extended)])
        self.assert_preflight_failed(result)
        self.assertIn("--allow-extended-ref", result.stderr)
        result = self.run_batch(extra=["--ref-hg38", str(extended), "--allow-extended-ref"])
        self.assertEqual(result.returncode, 0, result.stderr)

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

    def test_all_backends_share_dictionary_and_liftover_selection_and_options(self):
        self.install_fake_wrappers()
        for mode, name in (("picard-jar", "java"), ("gatk", "gatk"), ("picard-cli", "picard")):
            with self.subTest(mode=mode):
                self.out = self.root / f"result-{mode}"
                self.ref.with_suffix(".dict").unlink()
                trace = self.root / f"{mode}.jsonl"
                self.env["TEST_BACKEND_ARGS"] = str(trace)
                self.env.pop("JAVA_TOOL_OPTIONS", None)
                result = self.run_batch(extra=["--backend", mode, "--java-opts", "-Xmx64m -Dpattern=*"])
                self.assertEqual(result.returncode, 0, result.stderr)
                manifest = json.loads((self.out / "manifest.json").read_text())
                self.assertEqual(manifest["backend"]["kind"], mode)
                self.assertEqual(manifest["backend"]["java_options"], ["-Xmx64m", "-Dpattern=*"])
                calls = [json.loads(line) for line in trace.read_text().splitlines()]
                self.assertEqual(len(calls), 2)
                self.assertEqual({call["program"] for call in calls}, {name})
                self.assertIn("CreateSequenceDictionary", calls[0]["argv"])
                self.assertIn("LiftoverVcf", calls[1]["argv"])
                for call in calls:
                    if mode == "picard-cli":
                        self.assertEqual(call["JAVA_TOOL_OPTIONS"], '"-Xmx64m" "-Dpattern=*"')
                    elif mode == "gatk":
                        self.assertEqual(call["argv"][:2], ["--java-options", "-Xmx64m -Dpattern=*"])
                    else:
                        self.assertEqual(call["argv"][:2], ["-Xmx64m", "-Dpattern=*"])

    def test_explicit_jar_dictionary_does_not_use_available_gatk(self):
        self.install_fake_wrappers()
        self.ref.with_suffix(".dict").unlink()
        trace = self.root / "backend.jsonl"
        self.env["TEST_BACKEND_ARGS"] = str(trace)
        result = self.run_batch()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({json.loads(line)["program"] for line in trace.read_text().splitlines()}, {"java"})

    def test_jar_quoted_jvm_values_survive_dictionary_shell_bridge(self):
        self.ref.with_suffix(".dict").unlink()
        trace = self.root / "backend.jsonl"
        self.env["TEST_BACKEND_ARGS"] = str(trace)
        result = self.run_batch(extra=["--java-opts", '-Xmx64m -Dsynthetic="a b" -Dpattern=*'])
        self.assertEqual(result.returncode, 0, result.stderr)
        for line in trace.read_text().splitlines():
            self.assertEqual(json.loads(line)["argv"][:3], ["-Xmx64m", "-Dsynthetic=a b", "-Dpattern=*"])

    def test_dual_driver_uses_shared_backends_end_to_end(self):
        self.install_fake_wrappers()
        markers = self.root / "markers.vcf.gz"
        marker_header = ('##fileformat=VCFv4.2\n##contig=<ID=chrY,length=16>\n'
            '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">\n'
            '##INFO=<ID=HG,Number=1,Type=String,Description="Label">\n'
            '##INFO=<ID=ISOGG,Number=1,Type=String,Description="Label">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        with pysam.BGZFile(str(markers), "w") as handle:
            handle.write((marker_header + 'chrY\t2\tSYNTHETIC\tC\tG\t.\tPASS\tAA=C;HG=I-SYNTHETIC;ISOGG=I1\n').encode())
        pysam.tabix_index(str(markers), preset="vcf")
        gff = self.root / "markers.gff3"
        gff.write_text('##gff-version 3\nchrY\tpoint\tsnp\t2\t2\t.\t+\t.\t'
                       'ID=SYNTHETIC;allele_anc=C;allele_der=G;isogg_haplogroup=I1;yfull_node=I-SYNTHETIC\n')
        for mode in ("picard-jar", "gatk", "picard-cli"):
            with self.subTest(mode=mode):
                out = self.root / f"dual-{mode}"
                trace = self.root / f"dual-{mode}.jsonl"
                env = dict(self.env, TEST_MODE="success", TEST_BACKEND_ARGS=str(trace))
                result = subprocess.run(["sh", str(batch.SCRIPT_DIR / "run_y_dual_liftover_experiment.sh"),
                    "--sample-vcf-hs1", str(self.sample), "--ref-hs1", str(self.ref), "--ref-hg38", str(self.ref),
                    "--chain-hg38-to-hs1", str(self.chain), "--chain-hs1-to-hg38", str(self.chain),
                    "--markers-hg38", str(markers), "--markers-gff3-hg38", str(gff),
                    "--picard-jar", str(self.jar), "--backend", mode, "--out-dir", str(out)],
                    cwd=self.root, env=env, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue((out / "summary.tsv").is_file())
                for name in ("hs1_from_markers_lifted", "hg38_from_sample_lifted"):
                    rows = (out / f"{name}.marker_status.tsv").read_text().splitlines()
                    self.assertEqual(len(rows), 1)
                    self.assertEqual(rows[0].split("\t")[12], "derived")
                calls = [json.loads(line) for line in trace.read_text().splitlines()]
                self.assertEqual(len(calls), 2)
                expected = {"picard-jar": "java", "gatk": "gatk", "picard-cli": "picard"}[mode]
                self.assertEqual({call["program"] for call in calls}, {expected})

    def test_gatk_unrepresentable_jvm_value_fails_before_writes(self):
        self.install_fake_wrappers()
        Path(str(self.ref) + ".fai").unlink()
        result = self.run_batch(extra=["--backend", "gatk", "--java-opts", '-Xmx64m -Dsynthetic="a b"'])
        self.assert_preflight_failed(result)
        self.assertIn("GATK's launcher splits", result.stderr)
        self.assertFalse(Path(str(self.ref) + ".fai").exists())

    def test_dictionary_failure_keeps_stage_without_fallback_or_publication(self):
        self.install_fake_wrappers()
        self.ref.with_suffix(".dict").unlink()
        trace = self.root / "backend.jsonl"
        self.env["TEST_BACKEND_ARGS"] = str(trace)
        result = self.run_batch(mode="fail_dict")
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.out.exists())
        manifest = json.loads((self.stages()[0] / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual([json.loads(line)["program"] for line in trace.read_text().splitlines()], ["java"])

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

    def test_signals_stop_each_backend_tool_group_and_retain_stage(self):
        self.install_fake_wrappers()
        for mode, signum in ((m, s) for m in ("picard-jar", "gatk", "picard-cli")
                             for s in (signal.SIGTERM, signal.SIGINT)):
            with self.subTest(backend=mode, signal=signum):
                self.ready.unlink(missing_ok=True)
                proc = subprocess.Popen(self.command(extra=["--backend", mode]), env=dict(self.env, TEST_MODE="sleep"),
                                        cwd=self.root, text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
