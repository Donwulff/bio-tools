#!/usr/bin/env python3
"""Tiny synthetic FASTA/sidecar checks; no production reference downloads or writes."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference_utils import ReferenceError, dictionary_paths, read_fai

HERE = Path(__file__).resolve().parent
HAS_PYSAM = importlib.util.find_spec("pysam") is not None


@unittest.skipUnless(HAS_PYSAM, "indexed reference tests require pysam")
class ReferenceTest(unittest.TestCase):
    def setUp(self):
        import pysam
        self.tmp = tempfile.TemporaryDirectory(prefix="reference-check-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ref = self.root / "tiny.fa"
        self.ref.write_text(">chrY\nACGTACGT\n>chrM\nGGTT\n")
        pysam.faidx(str(self.ref))
        self.dictionary = self.root / "tiny.dict"
        self.write_dict()

    def write_dict(self):
        self.dictionary.write_text("@HD\tVN:1.6\n" + "".join(
            f"@SQ\tSN:{name}\tLN:{len(seq)}\tM5:{hashlib.md5(seq.encode()).hexdigest()}\n"
            for name, seq in [("chrY", "ACGTACGT"), ("chrM", "GGTT")]))

    def check(self, *args, ok=True):
        run = subprocess.run([sys.executable, str(HERE / "check_reference.py"), str(self.ref), *args],
                             text=True, capture_output=True)
        self.assertEqual(run.returncode == 0, ok, run.stderr)
        return run

    def test_valid_reference_and_md5(self):
        self.check()
        self.check("--check-md5")

    def test_dictionary_names_for_all_supported_suffixes(self):
        for suffix in (".fa", ".fna", ".fas", ".fasta", ".fa.gz", ".fna.gz", ".fas.gz", ".fasta.gz"):
            with self.subTest(suffix=suffix):
                paths = dictionary_paths(self.root / f"tiny{suffix}")
                self.assertEqual(paths[0], self.root / "tiny.dict")
                self.assertEqual(len(paths), 2 if suffix.endswith(".gz") else 1)

    def test_empty_fai_and_dictionary_are_rejected(self):
        for path in (Path(f"{self.ref}.fai"), self.dictionary):
            original = path.read_bytes()
            path.write_bytes(b"")
            self.check(ok=False)
            self.assertEqual(path.read_bytes(), b"")
            path.write_bytes(original)

    def test_dictionary_without_sq_is_rejected(self):
        self.dictionary.write_text("@HD\tVN:1.6\n")
        self.assertIn("@SQ", self.check(ok=False).stderr)

    def test_dictionary_names_lengths_order_and_duplicates_are_checked(self):
        original = self.dictionary.read_text()
        invalid = [original.replace("SN:chrY", "SN:other"),
                   original.replace("LN:8", "LN:9"),
                   "\n".join(original.splitlines()[::-1]) + "\n",
                   original + original.splitlines()[1] + "\n"]
        for text in invalid:
            with self.subTest(text=text):
                self.dictionary.write_text(text)
                self.check(ok=False)

    def test_nonempty_bad_fai_not_rebuilt(self):
        index = Path(f"{self.ref}.fai")
        index.write_text("invalid\n")
        before = index.read_bytes()
        self.check(ok=False)
        self.assertEqual(index.read_bytes(), before)

    def test_stale_plain_fasta_line_offsets_are_rejected(self):
        self.ref.write_text(self.ref.read_text().replace("ACGTACGT\n", "ACGT\nACGT\n"))
        self.check(ok=False)

    def test_invalid_and_duplicate_index_entries(self):
        path = Path(f"{self.ref}.fai")
        original = path.read_text()
        for text in (original + original.splitlines()[0] + "\n",
                     original.replace("\t8\t9\n", "\t8\t7\n"),
                     original.replace("chrY\t8", "chrY\t0")):
            with self.subTest(text=text):
                path.write_text(text)
                with self.assertRaises(ReferenceError):
                    read_fai(path)

    def test_same_length_sequence_change_needs_full_md5_check(self):
        self.ref.write_text(self.ref.read_text().replace("ACGTACGT", "ACGAACGT"))
        self.check()
        self.assertIn("M5", self.check("--check-md5", ok=False).stderr)

    def test_full_scan_requires_m5(self):
        self.dictionary.write_text("@SQ\tSN:chrY\tLN:8\n@SQ\tSN:chrM\tLN:4\n")
        self.check()
        self.assertIn("M5", self.check("--check-md5", ok=False).stderr)

    def test_plain_gzip_rejected_despite_gz_extension_and_sidecars(self):
        compressed = self.root / "tiny.fa.gz"
        with gzip.open(compressed, "wb") as fh:
            fh.write(self.ref.read_bytes())
        self.ref = compressed
        self.assertIn("not BGZF", self.check("--format-only", ok=False).stderr)
        self.assertFalse(Path(f"{self.ref}.fai").exists())

    def compress_reference(self):
        import pysam
        compressed = self.root / "tiny.fa.gz"
        pysam.tabix_compress(str(self.ref), str(compressed))
        pysam.faidx(str(compressed))
        self.ref = compressed

    def test_bgzf_and_aliases(self):
        self.compress_reference()
        self.check("--check-md5")
        alias = self.root / "tiny.fa.dict"
        alias.symlink_to(self.dictionary.name)
        self.check()
        alias.unlink()
        alias.write_text(self.dictionary.read_text().replace("LN:8", "LN:9"))
        self.check(ok=False)

    def test_conflicting_alias_checksums_are_rejected(self):
        self.compress_reference()
        alias = self.root / "tiny.fa.dict"
        digest = hashlib.md5(b"ACGTACGT").hexdigest()
        alias.write_text(self.dictionary.read_text().replace(digest, "0" * 32))
        self.assertIn("conflicting M5", self.check(ok=False).stderr)

    def test_broken_alias_is_not_ignored(self):
        self.compress_reference()
        (self.root / "tiny.fa.dict").symlink_to("missing.dict")
        self.check(ok=False)

    def test_gzi_presence_and_structure(self):
        self.compress_reference()
        index = Path(f"{self.ref}.gzi")
        original = index.read_bytes()
        for bad in (b"", b"bad", struct.pack("<Q", 10),
                    struct.pack("<QQQQQ", 2, 10, 20, 9, 30)):
            with self.subTest(bad=bad):
                index.write_bytes(bad)
                self.check(ok=False)
                self.assertEqual(index.read_bytes(), bad)
        index.write_bytes(original)
        self.check()

    @unittest.skipUnless(shutil.which("samtools"), "shell preparation integration requires samtools")
    def test_shared_shell_helper_creates_missing_sidecars_and_reuses_valid_ones(self):
        index = Path(f"{self.ref}.fai")
        index.unlink()
        self.dictionary.unlink()
        command = '. "$1/lib_liftover.sh"; SCRIPT_DIR="$1"; ensure_ref_sidecars "$2" "" ""'
        env = self.dictionary_tool_env()
        run = subprocess.run(["sh", "-eu", "-c", command, "test", str(HERE), str(self.ref)],
                             env=env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in (index, self.dictionary)}
        run = subprocess.run(["sh", "-eu", "-c", command, "test", str(HERE), str(self.ref)],
                             env=env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        for path, state in before.items():
            self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), state)

    def test_shell_helper_refuses_empty_dictionary_without_replacing_it(self):
        self.dictionary.write_bytes(b"")
        command = '. "$1/lib_liftover.sh"; SCRIPT_DIR="$1"; ensure_ref_sidecars "$2" "" ""'
        run = subprocess.run(["sh", "-eu", "-c", command, "test", str(HERE), str(self.ref)],
                             capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertEqual(self.dictionary.read_bytes(), b"")

    def dictionary_tool_env(self):
        tools = self.root / "bin"
        tools.mkdir(exist_ok=True)
        gatk = tools / "gatk"
        gatk.write_text('#!/bin/sh\n[ "$1" != --java-options ] || shift 2\nshift\nwhile [ "$#" -gt 0 ]; do\n'
                        'case "$1" in --REFERENCE) ref="$2"; shift 2;; --OUTPUT) out="$2"; shift 2;; *) exit 9;; esac\ndone\n'
                        'exec "$TEST_SAMTOOLS" dict -o "$out" "$ref"\n')
        gatk.chmod(0o755)
        env = {**os.environ, "PATH": f"{tools}:{os.environ['PATH']}", "TEST_SAMTOOLS": shutil.which("samtools")}
        for key in ("BIO_TOOLS_CFG", "Y_RESOURCES_DIR", "PICARD_JAR", "PICARD_JAVA_OPTS", "LIFTOVER_BACKEND",
                    "YREF_OUT_DIR", "YREF_HS1_NAME", "YREF_HG38_NAME", "YREF_PRIMARY_CONTIGS"):
            env.pop(key, None)
        return env

    def staging(self, *args, ok=True):
        run = subprocess.run(["sh", str(HERE / "prepare_y_refs.sh"),
                              "--ref-hs1-src", self.ref.name, "--ref-hg38-src", self.ref.name,
                              "--out-dir", "stage", "--hs1-name", "hs1.fa", "--hg38-name", "hg38.fa", *args],
                             cwd=self.root, capture_output=True, text=True, env=self.dictionary_tool_env())
        self.assertEqual(run.returncode == 0, ok, run.stderr)
        return run

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_invalid_backend_options_fail_before_staging(self):
        result = self.staging("--java-opts", '"unterminated', ok=False)
        self.assertIn("No closing quotation", result.stderr)
        self.assertFalse((self.root / "stage").exists())

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_unrepresentable_gatk_options_fail_before_copy_or_index_writes(self):
        self.dictionary.unlink()
        index = Path(str(self.ref) + ".fai")
        before = (index.read_bytes(), index.stat().st_mtime_ns)
        for mode in ("auto", "gatk"):
            with self.subTest(mode=mode):
                result = self.staging("--copy", "--backend", mode,
                                      "--java-opts", '-Dsynthetic="a b"', ok=False)
                self.assertIn("GATK's launcher splits", result.stderr)
                self.assertFalse((self.root / "stage").exists())
                self.assertEqual((index.read_bytes(), index.stat().st_mtime_ns), before)

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_noncopyable_source_dictionary_fails_before_staging_or_index_writes(self):
        self.dictionary.unlink()
        index = Path(str(self.ref) + ".fai")
        index.unlink()
        for kind in ("dangling", "directory"):
            with self.subTest(kind=kind):
                if kind == "dangling":
                    self.dictionary.symlink_to("missing.dict")
                else:
                    self.dictionary.mkdir()
                result = self.staging("--copy", "--backend", "gatk",
                                      "--java-opts", '-Dsynthetic="a b"', ok=False)
                self.assertIn("invalid dictionary path", result.stderr)
                self.assertFalse((self.root / "stage").exists())
                self.assertFalse(index.exists())
                if kind == "dangling":
                    self.assertTrue(self.dictionary.is_symlink())
                    self.dictionary.unlink()
                else:
                    self.assertTrue(self.dictionary.is_dir())
                    self.dictionary.rmdir()

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_relative_source_symlinks_are_absolute_and_reusable(self):
        self.staging()
        dest = self.root / "stage" / "hs1.fa"
        self.assertTrue(dest.is_symlink())
        self.assertEqual(dest.resolve(), self.ref)
        self.staging()

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_staging_copy_never_replaces_existing_target(self):
        self.staging("--copy")
        dest = self.root / "stage" / "hs1.fa"
        before = dest.read_bytes()
        self.staging("--copy", ok=False)
        self.assertEqual(dest.read_bytes(), before)

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_bgzf_staging_preserves_gzi_and_canonical_dictionary(self):
        self.compress_reference()
        self.staging("--copy")
        staged = self.root / "stage" / "hs1.fa"
        self.assertEqual(Path(f"{staged}.gzi").read_bytes(), Path(f"{self.ref}.gzi").read_bytes())
        self.assertEqual(staged.with_suffix(".dict").read_bytes(), self.dictionary.read_bytes())

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_orphaned_output_sidecar_is_not_reused(self):
        dest = self.root / "stage"
        dest.mkdir()
        (dest / "hs1.fa.fai").write_text("stale\n")
        self.staging(ok=False)
        self.assertFalse((dest / "hs1.fa").exists())

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_output_names_must_be_distinct(self):
        self.staging("--hg38-name", "hs1.fa", ok=False)
        self.staging("--hg38-name", "./hs1.fa", ok=False)
        self.assertFalse((self.root / "stage" / "hs1.fa").exists())

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_primary_subset_gets_its_own_index_layout(self):
        import pysam
        self.ref.write_text(">decoy\nAAAA\n>chrY\nACGTACGT\n>chrM\nGGTT\n")
        pysam.faidx(str(self.ref))
        self.staging("--primary-only", "--primary-contigs", "chrY chrM")
        staged = self.root / "stage" / "hs1.fa"
        self.assertNotIn(">decoy", staged.read_text())
        source_index = read_fai(Path(f"{self.ref}.fai"))
        target_index = read_fai(Path(f"{staged}.fai"))
        self.assertEqual(list(target_index), ["chrY", "chrM"])
        self.assertNotEqual(source_index["chrY"].offset, target_index["chrY"].offset)

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_failed_primary_extraction_leaves_no_partial_reference(self):
        self.dictionary_tool_env()
        fake = self.root / "bin" / "samtools"
        fake.write_text("#!/bin/sh\nexit 99\n")
        fake.chmod(0o755)
        self.staging("--primary-only", ok=False)
        self.assertFalse((self.root / "stage" / "hs1.fa").exists())
        self.assertFalse(list((self.root / "stage").glob("*.tmp.*")))

    def test_multiblock_bgzf_indexes_and_full_md5_scan(self):
        import pysam
        seq = "ACGT" * 20000
        self.ref.write_text(f">chrY\n{seq}\n>chrM\nGGTT\n")
        self.dictionary.write_text(
            f"@SQ\tSN:chrY\tLN:{len(seq)}\tM5:{hashlib.md5(seq.encode()).hexdigest()}\n"
            f"@SQ\tSN:chrM\tLN:4\tM5:{hashlib.md5(b'GGTT').hexdigest()}\n")
        pysam.faidx(str(self.ref))
        self.compress_reference()
        self.assertGreater(Path(f"{self.ref}.gzi").stat().st_size, 8)
        self.check("--check-md5")

    def test_crlf_and_last_line_without_newline(self):
        import pysam
        self.ref.write_bytes(b">chrY\r\nACGTACGT\r\n>chrM\r\nGGTT")
        pysam.faidx(str(self.ref))
        self.check("--check-md5")

    @unittest.skipUnless(shutil.which("samtools"), "staging integration requires samtools")
    def test_staging_handles_spaces_in_source_and_destination_paths(self):
        original = self.ref
        self.ref = self.root / "tiny reference.fa"
        original.rename(self.ref)
        Path(f"{original}.fai").rename(Path(f"{self.ref}.fai"))
        self.dictionary.rename(self.root / "tiny reference.dict")
        self.staging("--out-dir", "stage with spaces")
        self.assertEqual((self.root / "stage with spaces" / "hs1.fa").resolve(), self.ref)


if __name__ == "__main__":
    unittest.main()
