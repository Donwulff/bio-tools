#!/usr/bin/env python3
"""Backend selection/argv/Java transport tests with synthetic executables only."""

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

from liftover_backend import Backend, BackendError, java_options, select_backend, tool_options_value


HERE = Path(__file__).resolve().parent
FAKE_PROGRAM = '''import json, os, sys
print(json.dumps({"argv": sys.argv[1:], "JAVA_TOOL_OPTIONS": os.environ.get("JAVA_TOOL_OPTIONS")}))
sys.exit(int(os.environ.get("TEST_EXIT", "0")))
'''


class BackendTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="liftover-backend-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bin = self.root / "bin with spaces"
        self.bin.mkdir()
        self.repo = self.root / "repo"
        self.env = {"PATH": str(self.bin)}
        self.jar = self.root / "Picard with spaces.jar"
        self.jar.write_text("synthetic only")

    def exe(self, name):
        path = self.bin / name
        path.write_text(f"#!{sys.executable}\n" + FAKE_PROGRAM)
        path.chmod(0o755)
        return path

    def select(self, mode="auto", jar=None, options="-Xmx64m"):
        return select_backend(mode, jar, options, environ=self.env, repo_root=self.repo)

    def discovered_jar(self):
        path = self.repo / "mapping/picard.jar"
        path.parent.mkdir(parents=True)
        path.write_text("synthetic only")
        return path

    def test_explicit_jar_wins_auto_even_with_gatk(self):
        self.exe("java")
        self.exe("gatk")
        backend = self.select(jar=self.jar)
        self.assertEqual(backend.kind, "picard-jar")
        self.assertEqual(backend.jar, str(self.jar))

    def test_auto_gatk_precedes_discovered_jar_and_picard_cli(self):
        self.exe("gatk")
        self.exe("picard")
        self.discovered_jar()
        self.assertEqual(self.select().kind, "gatk")

    def test_auto_discovered_jar_precedes_picard_cli(self):
        self.exe("java")
        self.exe("picard")
        path = self.discovered_jar()
        backend = self.select()
        self.assertEqual(backend.kind, "picard-jar")
        self.assertEqual(backend.jar, str(path))

    def test_auto_picard_cli_fallback(self):
        self.exe("picard")
        self.assertEqual(self.select().kind, "picard-cli")

    def test_missing_or_empty_supplied_jar_never_falls_back(self):
        self.exe("gatk")
        self.exe("java")
        empty = self.root / "empty.jar"
        empty.touch()
        for path in (empty, self.root / "missing.jar"):
            with self.subTest(path=path), self.assertRaisesRegex(BackendError, "Picard jar"):
                self.select(jar=path)

    def test_explicit_nonjar_mode_ignores_global_jar(self):
        for mode, executable in (("gatk", "gatk"), ("picard-cli", "picard")):
            with self.subTest(mode=mode):
                self.exe(executable)
                self.assertIsNone(self.select(mode, jar=self.root / "missing.jar").jar)

    def test_explicit_jar_mode_requires_jar_and_java(self):
        self.exe("java")
        with self.assertRaisesRegex(BackendError, "no Picard jar"):
            self.select("picard-jar")
        (self.bin / "java").unlink()
        self.exe("gatk")
        with self.assertRaisesRegex(BackendError, "missing backend command: java"):
            self.select("picard-jar", self.jar)

    def test_missing_requested_wrapper_does_not_fall_back(self):
        self.exe("java")
        with self.assertRaisesRegex(BackendError, "missing backend command: gatk"):
            self.select("gatk", self.jar)

    def test_unknown_backend_or_tool_rejected(self):
        with self.assertRaisesRegex(BackendError, "unknown backend"):
            self.select("unknown")
        self.exe("java")
        with self.assertRaisesRegex(BackendError, "unsupported Picard tool"):
            self.select(jar=self.jar).invocation("OtherTool", [])

    def test_pinned_backend_does_not_reselect_after_path_change(self):
        self.exe("picard")
        selected = self.select()
        self.exe("gatk")
        pinned = Backend.from_metadata(json.loads(json.dumps(selected.metadata())))
        self.assertEqual(pinned, selected)
        self.assertEqual(self.select().kind, "gatk")

    def test_invalid_pinned_metadata_is_rejected(self):
        self.exe("java")
        data = self.select(jar=self.jar).metadata()
        for bad in ({}, {**data, "kind": "auto"}, {**data, "executable": "relative"},
                    {**data, "java_options": "-Xmx64m"}, {**data, "jar": None},
                    {**data, "java_options": [7]}, {**data, "extra": True},
                    {**data, "kind": "gatk"}):
            with self.subTest(bad=bad), self.assertRaises(BackendError):
                Backend.from_metadata(bad)

    def test_missing_pinned_executable_fails(self):
        executable = self.exe("java")
        data = self.select(jar=self.jar).metadata()
        executable.unlink()
        with self.assertRaisesRegex(BackendError, "missing/not executable"):
            Backend.from_metadata(data)

    def test_jvm_parser_keeps_quoted_values_and_literals(self):
        self.assertEqual(java_options('-Xmx64m -Dsynthetic="path with spaces" -Dpattern=*'),
                         ("-Xmx64m", "-Dsynthetic=path with spaces", "-Dpattern=*"))

    def test_malformed_jvm_options_rejected(self):
        for raw in ('"unterminated', "className", "-Xmx64m 64m", "-jar -version", "--module=other/main",
                    "--class-path=other.jar", "-Dk='x\ny'", "-Dk=x\0"):
            with self.subTest(raw=raw), self.assertRaises(BackendError):
                java_options(raw)

    def test_explicit_empty_jvm_options_preserved(self):
        self.exe("java")
        self.assertEqual(self.select(jar=self.jar, options="").java_options, ())

    def test_jar_argv_preserves_spaces_and_never_expands_globs(self):
        self.exe("java")
        backend = self.select(jar=self.jar, options='-Dsynthetic="path with spaces" -Dpattern=*')
        arguments = ["--INPUT", "input with spaces.vcf", "--OUTPUT", "literal$(touch BAD)*.vcf"]
        argv, updates = backend.invocation("LiftoverVcf", arguments)
        self.assertEqual(updates, {})
        self.assertEqual(argv[1:3], ["-Dsynthetic=path with spaces", "-Dpattern=*"])
        self.assertEqual(argv[-len(arguments):], arguments)
        run = subprocess.run(argv, env=self.env, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(json.loads(run.stdout)["argv"], argv[1:])
        self.assertFalse((self.root / "BAD").exists())

    def test_gatk_jvm_option_transport_and_tool_paths(self):
        self.exe("gatk")
        argv, updates = self.select("gatk", options="-Xmx64m -Dpattern=*").invocation(
            "CreateSequenceDictionary", ["--REFERENCE", "reference with spaces.fa"])
        self.assertEqual(argv[1:], ["--java-options", "-Xmx64m -Dpattern=*",
                                   "CreateSequenceDictionary", "--REFERENCE", "reference with spaces.fa"])
        self.assertEqual(updates, {})

    def test_gatk_rejects_space_containing_jvm_value_instead_of_splitting(self):
        self.exe("gatk")
        with self.assertRaisesRegex(BackendError, "GATK's launcher splits"):
            self.select("gatk", options='-Dsynthetic="path with spaces"')

    def test_picard_cli_options_append_without_mutating_environment(self):
        self.exe("picard")
        env = {"JAVA_TOOL_OPTIONS": "-Dexisting=kept"}
        argv, updates = self.select("picard-cli", options='-Xmx64m -Dsynthetic="a b"').invocation(
            "LiftoverVcf", ["--INPUT", "a b.vcf"], env)
        self.assertEqual(argv[1:], ["LiftoverVcf", "--INPUT", "a b.vcf"])
        self.assertEqual(updates, {"JAVA_TOOL_OPTIONS": '-Dexisting=kept "-Xmx64m" "-Dsynthetic=a b"'})
        self.assertEqual(env, {"JAVA_TOOL_OPTIONS": "-Dexisting=kept"})

    def test_picard_cli_empty_options_do_not_clear_existing_options(self):
        self.exe("picard")
        _argv, updates = self.select("picard-cli", options="").invocation("LiftoverVcf", [])
        self.assertEqual(updates, {})

    def test_unrepresentable_mixed_quotes_rejected_for_picard_cli(self):
        self.exe("picard")
        with self.assertRaisesRegex(BackendError, "mixed quote"):
            self.select("picard-cli", options=shlex.join(['-Dsynthetic=a\'"b']))

    @unittest.skipUnless(shutil.which("java"), "optional runtime option transport check requires Java")
    def test_java_runtime_preserves_environment_option_quotes(self):
        for value in ("a b", "a'b c", 'a"b c', r"path\with\backslashes", "literal$(noexec)*"):
            with self.subTest(value=value):
                env = {k: v for k, v in os.environ.items()
                       if k not in ("JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS")}
                env["JAVA_TOOL_OPTIONS"] = tool_options_value(("-Xmx64m", "-Dsynthetic=" + value))
                run = subprocess.run([shutil.which("java"), "-XshowSettings:properties", "-version"],
                                     env=env, capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
                self.assertIn("    synthetic = " + value + "\n", run.stderr)

    def test_describe_selects_without_executing_backend(self):
        self.exe("java")
        command = [sys.executable, str(HERE / "liftover_backend.py"), "--describe",
                   "--picard-jar", str(self.jar), "--java-opts=-Xmx64m"]
        run = subprocess.run(command, env=self.env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        data = json.loads(run.stdout)
        self.assertEqual(data["kind"], "picard-jar")
        self.assertNotIn("argv", data)

    def test_option_syntax_check_requires_no_backend(self):
        run = subprocess.run([sys.executable, str(HERE / "liftover_backend.py"),
                              "--check-options", "--backend", "gatk", "--java-opts=-Xmx64m"],
                             env=self.env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(run.stdout, "")

    def test_cli_exec_preserves_backend_failure_and_arguments(self):
        self.exe("picard")
        self.env.update(LIFTOVER_BACKEND="picard-cli", PICARD_JAVA_OPTS="-Xmx64m", TEST_EXIT="23")
        run = subprocess.run([sys.executable, str(HERE / "liftover_backend.py"), "LiftoverVcf",
                              "--INPUT", "a b.vcf"], env=self.env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 23, run.stderr)
        data = json.loads(run.stdout)
        self.assertEqual(data["argv"], ["LiftoverVcf", "--INPUT", "a b.vcf"])
        self.assertEqual(data["JAVA_TOOL_OPTIONS"], '"-Xmx64m"')

    def test_shell_helper_runs_pinned_backend_for_all_modes(self):
        # Shell helpers need only sh + Python here; the tool executables are fakes.
        for name in ("java", "gatk", "picard"):
            self.exe(name)
        self.env["PATH"] += os.pathsep + os.environ["PATH"]
        command = ('. "$1/lib_liftover.sh"; SCRIPT_DIR="$1"; '
                   'spec=$(select_liftover_backend "$2" "$3" "-Xmx64m"); '
                   'run_liftover_tool "$spec" CreateSequenceDictionary --REFERENCE "a b.fa"')
        for mode in ("gatk", "picard-jar", "picard-cli"):
            with self.subTest(mode=mode):
                run = subprocess.run(["sh", "-eu", "-c", command, "test", str(HERE), mode, str(self.jar)],
                                     env=self.env, capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
                data = json.loads(run.stdout)
                self.assertEqual(data["argv"][-3:], ["CreateSequenceDictionary", "--REFERENCE", "a b.fa"])


if __name__ == "__main__":
    unittest.main()
