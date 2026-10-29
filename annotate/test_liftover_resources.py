#!/usr/bin/env python3
"""Synthetic discovery/configuration checks; no network or private inputs."""

import contextlib
import io
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

from liftover_resources import (
    EXTENDED_NAMES, REFERENCE_NAMES, ResourceError, config_environment,
    find_reference, reference_roots,
)


HERE = Path(__file__).resolve().parent
CONFIG_KEYS = ("BIO_TOOLS_CFG", "REF_ROOT", "REF_HS1", "REF_HG38", "BIO_TOOLS_INDEX_DIR",
               "INDEX_DIR", "Y_RESOURCES_DIR", "YMARKERS_OUT", "YMARKERS_GFF3_OUT",
               "PICARD_JAR", "PICARD_JAVA_OPTS", "YREF_OUT_DIR", "YREF_HS1_NAME", "YREF_HG38_NAME")


class ResourceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="liftover-resources-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.cwd = self.root / "analysis with spaces"
        self.repo = self.root / "source repo"
        self.cwd.mkdir()
        self.repo.mkdir()
        self.env = {k: v for k, v in os.environ.items() if k not in CONFIG_KEYS}

    def fixture(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(">chrY\nACGT\n")
        return path

    def find(self, assembly="hg38", **kwargs):
        return find_reference(assembly, environ=self.env, cwd=self.cwd, repo_root=self.repo, **kwargs)

    def config(self, text, name="bio-tools.cfg"):
        path = self.cwd / name
        path.write_text(text)
        return path

    def test_explicit_reference_wins_over_environment_and_search(self):
        chosen = self.fixture(self.cwd / "chosen.fa")
        self.env["REF_HG38"] = "missing.fa"
        self.fixture(self.cwd / "resources/ref" / REFERENCE_NAMES["hg38"][0])
        self.assertEqual(self.find(explicit=chosen), chosen)

    def test_missing_or_empty_explicit_reference_does_not_fall_back(self):
        self.fixture(self.cwd / "resources/ref" / REFERENCE_NAMES["hg38"][0])
        for path in (self.cwd / "missing.fa", self.cwd / "empty.fa"):
            with self.subTest(path=path):
                if path.name == "empty.fa":
                    path.touch()
                with self.assertRaisesRegex(ResourceError, "readable, nonempty"):
                    self.find(explicit=path)

    def test_reference_environment_values_work_for_both_assemblies(self):
        for assembly in REFERENCE_NAMES:
            path = self.fixture(self.cwd / f"{assembly} chosen.fa")
            self.env[f"REF_{assembly.upper()}"] = path.name
            self.assertEqual(self.find(assembly), path)

    def test_root_precedence_and_spaces(self):
        roots = {key: self.root / key for key in ("REF_ROOT", "BIO_TOOLS_INDEX_DIR", "INDEX_DIR")}
        for key, root in roots.items():
            self.env[key] = str(root)
            self.fixture(root / REFERENCE_NAMES["hg38"][0])
        chosen = self.fixture(self.root / "CLI root with spaces" / REFERENCE_NAMES["hg38"][0])
        self.assertEqual(self.find(root=chosen.parent), chosen)
        for key, root in roots.items():
            self.assertEqual(self.find(), root / REFERENCE_NAMES["hg38"][0])
            del self.env[key]

    def test_relative_configured_root_is_analysis_relative(self):
        chosen = self.fixture(self.cwd / "relative root" / REFERENCE_NAMES["hg38"][0])
        self.assertEqual(self.find(root="relative root"), chosen)

    def test_resource_root_before_working_directory_and_repo_fallbacks(self):
        external = self.root / "shared resources"
        self.env["Y_RESOURCES_DIR"] = str(external)
        roots = [external / "ref", self.cwd / "resources/ref", self.repo / "mapping/index"]
        paths = [self.fixture(root / REFERENCE_NAMES["hg38"][0]) for root in roots]
        for path in paths:
            self.assertEqual(self.find(), path)
            path.unlink()

    def test_roots_are_not_split_or_duplicated(self):
        root = self.cwd / "resources/ref"
        self.env["Y_RESOURCES_DIR"] = str(root.parent)
        self.assertEqual(reference_roots(str(root), environ=self.env, cwd=self.cwd, repo_root=self.repo),
                         [root, self.repo / "mapping/index"])

    def test_plain_fasta_precedes_compressed_within_a_root(self):
        for assembly in REFERENCE_NAMES:
            base = self.cwd / assembly
            plain = self.fixture(base / REFERENCE_NAMES[assembly][0])
            compressed = self.fixture(base / REFERENCE_NAMES[assembly][1])
            self.assertEqual(self.find(assembly, root=base), plain)
            plain.unlink()
            # Selection is not a claim that gzip content is BGZF or that sidecars are valid.
            self.assertEqual(self.find(assembly, root=base), compressed)

    def test_staging_generic_hg38_name_is_discoverable(self):
        chosen = self.fixture(self.cwd / "resources/ref/hg38.fa")
        self.assertEqual(self.find(), chosen)

    def test_extended_names_are_not_automatic_defaults(self):
        for assembly, names in EXTENDED_NAMES.items():
            for name in names:
                with self.subTest(assembly=assembly, name=name):
                    chosen = self.fixture(self.cwd / "resources/ref" / name)
                    with self.assertRaisesRegex(ResourceError, "auto-detect failed"):
                        self.find(assembly)
                    with self.assertRaisesRegex(ResourceError, "allow-extended-ref"):
                        self.find(assembly, explicit=chosen)
                    self.assertEqual(self.find(assembly, allow_extended=True), chosen)
                    chosen.unlink()

    def test_base_reference_wins_across_roots_even_if_extended_allowed(self):
        extended = self.fixture(self.root / "index" / EXTENDED_NAMES["hg38"][0])
        base = self.fixture(self.cwd / "resources/ref" / REFERENCE_NAMES["hg38"][0])
        self.assertEqual(self.find(root=extended.parent, allow_extended=True), base)
        self.assertEqual(self.find(explicit=extended, allow_extended=True), extended)

    def test_neutral_symlink_name_does_not_hide_known_extended_target(self):
        target = self.fixture(self.root / "custom.IMGT.fa")
        alias = self.cwd / "neutral.fa"
        alias.symlink_to(target)
        with self.assertRaisesRegex(ResourceError, "extended/decoy"):
            self.find(explicit=alias)
        self.assertEqual(self.find(explicit=alias, allow_extended=True), alias)

    def test_broken_discovered_symlink_is_not_ignored(self):
        path = self.cwd / "resources/ref" / REFERENCE_NAMES["hg38"][0]
        path.parent.mkdir(parents=True)
        path.symlink_to("missing.fa")
        with self.assertRaisesRegex(ResourceError, "readable, nonempty"):
            self.find()

    def test_no_recursive_filesystem_search(self):
        self.fixture(self.repo / "somewhere else" / REFERENCE_NAMES["hg38"][0])
        with self.assertRaisesRegex(ResourceError, "auto-detect failed"):
            self.find()

    def test_unknown_assembly_rejected(self):
        with self.assertRaisesRegex(ResourceError, "unsupported assembly"):
            self.find("other")

    def test_config_optional_default_and_required_explicit_file(self):
        self.assertEqual(config_environment(self.env, self.cwd)["PATH"], self.env["PATH"])
        self.env["BIO_TOOLS_CFG"] = "missing.cfg"
        with contextlib.redirect_stderr(io.StringIO()) as errors:
            with self.assertRaisesRegex(ResourceError, "configuration"):
                config_environment(self.env, self.cwd)
        self.assertIn("BIO_TOOLS_CFG is not a file", errors.getvalue())

    def test_config_overrides_environment_without_requiring_export(self):
        self.env["REF_ROOT"] = "old"
        self.config("REF_ROOT='shared index with spaces'\nY_RESOURCES_DIR=resources2\n")
        loaded = config_environment(self.env, self.cwd)
        self.assertEqual(loaded["REF_ROOT"], "shared index with spaces")
        self.assertEqual(loaded["Y_RESOURCES_DIR"], "resources2")

    def test_explicit_config_basename_is_sourced_from_analysis_directory(self):
        self.config("REF_HG38='custom ref.fa'\n", name="local.cfg")
        self.env["BIO_TOOLS_CFG"] = "local.cfg"
        self.assertEqual(config_environment(self.env, self.cwd)["REF_HG38"], "custom ref.fa")

    def test_config_output_is_diagnostic_not_environment_data(self):
        self.config("echo config-notice\nREF_ROOT='a=b\nc d'\n")
        with contextlib.redirect_stderr(io.StringIO()) as errors:
            loaded = config_environment(self.env, self.cwd)
        self.assertEqual(loaded["REF_ROOT"], "a=b\nc d")
        self.assertIn("config-notice", errors.getvalue())

    def test_config_failure_is_not_ignored(self):
        self.config("return 7\n")
        with self.assertRaisesRegex(ResourceError, "configuration"):
            config_environment(self.env, self.cwd)

    def test_config_can_unset_an_inherited_override(self):
        self.env["REF_ROOT"] = "old"
        self.config("unset REF_ROOT\n")
        self.assertNotIn("REF_ROOT", config_environment(self.env, self.cwd))

    def test_shell_loader_restores_autoexport_state(self):
        self.config("REF_ROOT=chosen\n")
        for enabled in (False, True):
            command = ('. "$1/lib_liftover.sh"; ' + ("set -a; " if enabled else "")
                       + 'load_liftover_config; printf "%s" "$-"')
            run = subprocess.run(["sh", "-eu", "-c", command, "test", str(HERE)],
                                 cwd=self.cwd, env=self.env, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual("a" in run.stdout, enabled)

    def test_shell_discovery_uses_python_selector_with_spaces(self):
        chosen = self.fixture(self.cwd / "shared index" / REFERENCE_NAMES["hs1"][0])
        command = '. "$1/lib_liftover.sh"; SCRIPT_DIR="$1"; find_liftover_ref hs1 "" "$2" 0'
        run = subprocess.run(["sh", "-eu", "-c", command, "test", str(HERE), str(chosen.parent)],
                             cwd=self.cwd, env=self.env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(run.stdout.strip(), str(chosen))

    @unittest.skipUnless(all(shutil.which(t) for t in ("java", "bcftools", "tabix", "samtools")),
                         "dual-driver smoke check requires installed commands and pysam")
    def test_dual_driver_cli_overrides_config_before_selection(self):
        try:
            import pysam  # noqa: F401
        except ImportError:
            self.skipTest("requires pysam")
        hs1 = self.fixture(self.cwd / "chosen hs1.fa")
        hg38 = self.fixture(self.cwd / "chosen hg38.fa")
        self.config("REF_HS1=missing.fa\nREF_HG38=missing.fa\n")
        chain = self.cwd / "fake.chain"
        chain.write_text("synthetic only")
        sample = self.cwd / "sample.vcf"
        sample.write_text("synthetic only")
        run = subprocess.run(["sh", str(HERE / "run_y_dual_liftover_experiment.sh"),
                              "--sample-vcf-hs1", str(sample), "--chain-hg38-to-hs1", str(chain),
                              "--chain-hs1-to-hg38", str(chain), "--ref-hs1", str(hs1),
                              "--ref-hg38", str(hg38), "--picard-jar", "missing.jar"],
                             cwd=self.cwd, env=self.env, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn(f"hs1  {hs1}", run.stdout)
        self.assertIn(f"hg38 {hg38}", run.stdout)
        self.assertIn("missing.jar", run.stderr)
        self.assertFalse(Path(str(hs1) + ".fai").exists())

    def test_chain_fetcher_uses_configured_root_or_cli_without_network(self):
        tools = self.root / "bin"
        tools.mkdir()
        wget = tools / "wget"
        wget.write_text('#!/bin/sh\nwhile [ "$#" -gt 0 ]; do\n'
                        'case "$1" in -O) out="$2"; shift 2;; *) shift;; esac\ndone\n'
                        'printf "synthetic chain\\n" > "$out"\n')
        wget.chmod(0o755)
        self.env["PATH"] = str(tools) + os.pathsep + self.env["PATH"]
        shared = self.root / "shared resources"
        self.config(f"Y_RESOURCES_DIR={shlex.quote(str(shared))}\n")
        for out, args in ((shared / "chains", []),
                          (self.cwd / "CLI chains", ["--out-dir", "CLI chains"])):
            with self.subTest(out=out):
                run = subprocess.run(["sh", str(HERE / "fetch_liftover_chains.sh"), *args],
                                     cwd=self.cwd, env=self.env, capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertEqual(len(list(out.glob("*.chain.gz"))), 2)
        self.assertFalse((self.cwd / "resources").exists())

    @unittest.skipUnless(all(shutil.which(t) for t in ("bcftools", "tabix", "bgzip")),
                         "marker fetcher smoke check requires htslib tools and pysam")
    def test_marker_fetcher_reuses_configured_or_cli_resources_without_network(self):
        try:
            import pysam
        except ImportError:
            self.skipTest("requires pysam")
        tools = self.root / "bin"
        tools.mkdir()
        wget = tools / "wget"
        wget.write_text('#!/bin/sh\necho "unexpected download" >&2\nexit 99\n')
        wget.chmod(0o755)
        self.env["PATH"] = str(tools) + os.pathsep + self.env["PATH"]
        shared = self.root / "shared resources"
        self.config(f"Y_RESOURCES_DIR={shlex.quote(str(shared))}\n")
        header = ('##fileformat=VCFv4.2\n##contig=<ID=chrY,length=4>\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for out, args in ((shared / "snps_hg38.vcf.gz", []),
                          (self.cwd / "CLI markers.vcf.gz", ["--out", "CLI markers.vcf.gz"])):
            with self.subTest(out=out):
                out.parent.mkdir(parents=True, exist_ok=True)
                with pysam.BGZFile(str(out), "w") as handle:
                    handle.write((header + 'chrY\t2\tSYNTHETIC\tC\tG\t.\tPASS\t.\n').encode())
                pysam.tabix_index(str(out), preset="vcf")
                out.with_name(out.name.removesuffix(".vcf.gz") + ".gff3").write_text("##gff-version 3\n")
                run = subprocess.run(["sh", str(HERE / "fetch_ybrowse_markers.sh"), *args],
                                     cwd=self.cwd, env=self.env, capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertNotIn("unexpected download", run.stderr)
                if not args:
                    self.assertIn(str(out), run.stdout)
        self.assertFalse((self.cwd / "resources").exists())

    @unittest.skipUnless(shutil.which("samtools"), "reference staging smoke check requires samtools and pysam")
    def test_staging_defaults_follow_resource_root_and_are_discoverable(self):
        try:
            import pysam
        except ImportError:
            self.skipTest("requires pysam")
        source = self.fixture(self.cwd / "synthetic.fa")
        pysam.faidx(str(source))
        source.with_suffix(".dict").write_text("@SQ\tSN:chrY\tLN:4\n")
        shared = self.root / "shared resources"
        self.config(f"Y_RESOURCES_DIR={shlex.quote(str(shared))}\n")
        run = subprocess.run(["sh", str(HERE / "prepare_y_refs.sh"),
                              "--ref-hs1-src", str(source), "--ref-hg38-src", str(source)],
                             cwd=self.cwd, env=self.env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        loaded = config_environment(self.env, self.cwd)
        for assembly, name in (("hs1", REFERENCE_NAMES["hs1"][0]), ("hg38", "hg38.fa")):
            chosen = find_reference(assembly, environ=loaded, cwd=self.cwd, repo_root=self.repo)
            self.assertEqual(chosen, shared / "ref" / name)
            self.assertEqual(chosen.resolve(), source)
        self.assertFalse((self.cwd / "resources").exists())


if __name__ == "__main__":
    unittest.main()
