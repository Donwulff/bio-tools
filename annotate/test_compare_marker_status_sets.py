#!/usr/bin/env python3
"""Synthetic comparator regressions. Run this file directly; no sample inputs."""
from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent


def row(pos=1, marker="SYNTHETIC1", status="derived", ref="A", alt="C", aa="A", reason="."):
    return ["chrY", str(pos), marker, ref, alt, aa, "G-SYNTHETIC", "G2",
            "1" if status == "derived" else "0", "40", "10", "variant", status, reason]


class ComparatorTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="marker-comparator-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.a, self.b = self.root / "a.tsv", self.root / "b.tsv"
        self.out = self.root / "comparison"
        self.write(self.a, [row()])
        self.write(self.b, [row()])

    @staticmethod
    def write(path, rows):
        with path.open("w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)

    def run_cli(self, *args, ok=True, inputs=None):
        cmd = [sys.executable, str(HERE / "compare_marker_status_sets.py")]
        for name, path in inputs or [("a", self.a), ("b", self.b)]:
            cmd.extend(["--input", f"{name}={path}"])
        cmd.extend(["--out-prefix", str(self.out), *args])
        result = subprocess.run(cmd, text=True, capture_output=True, cwd=self.root)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def report(self, suffix):
        with Path(f"{self.out}.{suffix}.tsv").open() as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def assert_no_reports(self):
        self.assertFalse(list(self.root.glob("comparison.*.tsv")))

    def test_conflicting_duplicate_rejected_before_output(self):
        self.write(self.a, [row(), row(status="ancestral")])
        result = self.run_cli(ok=False)
        self.assertIn("duplicate", result.stderr.lower())
        self.assert_no_reports()

    def test_identical_duplicate_requires_explicit_policy(self):
        self.write(self.a, [row(), row()])
        self.run_cli(ok=False)
        result = self.run_cli("--duplicate-policy", "collapse-identical")
        self.assertIn("collapsed", result.stderr.lower())
        self.assertEqual(self.report("pairwise")[0]["callable_both"], "1")

    def test_conflicting_duplicate_never_collapsed(self):
        self.write(self.a, [row(), row(status="ancestral")])
        self.run_cli("--duplicate-policy", "collapse-identical", ok=False)
        self.assert_no_reports()

    def test_duplicate_dataset_names_are_rejected(self):
        result = self.run_cli(inputs=[("same", self.a), ("same", self.b)], ok=False)
        self.assertIn("dataset", result.stderr.lower())
        self.assert_no_reports()

    def test_empty_dataset_path_is_rejected(self):
        self.run_cli(inputs=[("a", ""), ("b", self.b)], ok=False)

    def test_incompatible_alleles_rejected_even_when_statuses_agree(self):
        self.write(self.b, [row(alt="G")])
        result = self.run_cli(ok=False)
        self.assertIn("allele", result.stderr.lower())
        self.assert_no_reports()

    def test_ref_alt_swap_with_same_ancestral_derived_alleles_is_compatible(self):
        self.write(self.b, [row(ref="C", alt="A", aa="A")])
        self.run_cli()
        self.assertEqual(self.report("pairwise")[0]["callable_same"], "1")

    def test_complement_requires_explicit_cross_assembly_option(self):
        self.write(self.b, [row(pos=8, ref="T", alt="G", aa="T")])
        self.run_cli("--key-mode", "id", ok=False)
        self.run_cli("--key-mode", "id", "--allow-reverse-complement")
        self.assertEqual(self.report("pairwise")[0]["reverse_complemented_keys"], "1")

    def test_reverse_complement_option_not_allowed_for_same_locus_keys(self):
        self.run_cli("--allow-reverse-complement", ok=False)

    def test_same_id_at_multiple_positions_rejected_in_id_mode(self):
        self.write(self.a, [row(), row(pos=2)])
        self.run_cli("--key-mode", "id", ok=False)

    def test_anonymous_id_rejected_in_id_mode(self):
        self.write(self.a, [row(marker=".")])
        self.run_cli("--key-mode", "id", ok=False)

    def test_conflicting_alias_calls_at_same_locus_are_rejected(self):
        self.write(self.a, [row(), row(marker="SYNTHETIC_ALIAS", status="ancestral")])
        self.run_cli(ok=False)
        self.assert_no_reports()

    def test_aliases_do_not_inflate_unique_locus_count(self):
        rows = [row(), row(marker="SYNTHETIC_ALIAS")]
        self.write(self.a, rows)
        self.write(self.b, rows)
        self.run_cli()
        report = self.report("pairwise")[0]
        self.assertEqual(report["callable_both"], "2")
        self.assertEqual(report["callable_unique_locus_pairs"], "1")

    def test_callable_denominator_excludes_matching_nocalls(self):
        common = [row(i, f"SYNTHETIC{i}", "nocall") for i in range(1, 4)]
        self.write(self.a, common + [row(4, "SYNTHETIC4")])
        self.write(self.b, common + [row(4, "SYNTHETIC4", "ancestral")])
        self.run_cli()
        report = self.report("pairwise")[0]
        self.assertEqual(report["same_status"], "3")
        self.assertEqual(report["callable_both"], "1")
        self.assertEqual(report["callable_different"], "1")
        self.assertEqual(report["callable_agreement_pct"], "0.000000")

    def test_zero_callable_denominator_is_na(self):
        self.write(self.a, [row(status="filtered", reason="low_gq")])
        self.write(self.b, [row(status="nocall", reason="missing_gt")])
        self.run_cli()
        self.assertEqual(self.report("pairwise")[0]["callable_agreement_pct"], "NA")
        self.assertEqual(self.report("matrix")[0]["a.reason"], "low_gq")

    def test_third_dataset_missing_pairs_do_not_inflate_callability(self):
        third = self.root / "c.tsv"
        self.write(third, [row(2, "SYNTHETIC2")])
        self.run_cli(inputs=[("a", self.a), ("b", self.b), ("c", third)])
        report = self.report("pairwise")[0]
        self.assertEqual(report["keys_global_union"], "2")
        self.assertEqual(report["keys_pair_union"], "1")
        self.assertEqual(report["callable_both"], "1")

    def test_malformed_rows_are_not_silently_skipped(self):
        self.write(self.a, [["chrY", "1", "SYNTHETIC1"]])
        self.run_cli(ok=False)
        self.assert_no_reports()

    def test_invalid_positions_and_states_are_rejected(self):
        for bad in (row(pos=0), row(pos="bad"), row(status="typo")):
            with self.subTest(bad=bad):
                self.write(self.a, [bad])
                self.run_cli(ok=False)

    def test_comment_header_and_legacy_13_columns_are_supported(self):
        self.write(self.a, [["# comment"],
                           ["CHROM", "POS", "ID", "REF", "ALT", "AA", "HG", "ISOGG", "GT", "GQ", "DP", "SOURCE", "STATUS"],
                           row()[:13]])
        self.run_cli()
        self.assertEqual(self.report("pairwise")[0]["callable_same"], "1")

    def test_output_prefix_suffix_is_not_lost(self):
        self.out = self.root / "run.version1"
        self.run_cli()
        self.assertTrue(Path(f"{self.out}.matrix.tsv").exists())
        self.assertFalse((self.root / "run.matrix.tsv").exists())

    def test_existing_output_refused_unless_force(self):
        self.run_cli()
        self.run_cli(ok=False)
        self.run_cli("--force")

    def test_expected_panel_includes_missing_and_filtered_markers(self):
        expected = self.root / "expected.tsv"
        expected.write_text("id\tstatus\nSYNTHETIC1\tderived\nSYNTHETIC2\tancestral\nSYNTHETIC3\tderived\n")
        self.write(self.a, [row(), row(2, "SYNTHETIC2", "filtered", reason="low_dp")])
        self.write(self.b, [row(status="ancestral")])
        self.run_cli("--expected-states", str(expected))
        details = self.report("expected")
        a = {r["ID"]: r for r in details if r["dataset"] == "a"}
        b = {r["ID"]: r for r in details if r["dataset"] == "b"}
        self.assertEqual(a["SYNTHETIC1"]["comparison"], "match")
        self.assertEqual(a["SYNTHETIC2"]["status"], "filtered")
        self.assertEqual(a["SYNTHETIC3"]["status"], "missing")
        self.assertEqual(b["SYNTHETIC1"]["comparison"], "conflict")

    def test_expected_panel_conflicts_and_excluded_ids_are_rejected(self):
        expected = self.root / "expected.tsv"
        expected.write_text("SYNTHETIC1\tderived\nSYNTHETIC1\tancestral\n")
        self.run_cli("--expected-states", str(expected), ok=False)
        expected.write_text("SYNTHETIC1\tderived\n")
        self.run_cli("--expected-states", str(expected), "--id-regex", "^NOT_INCLUDED$", ok=False)

    def test_expected_id_at_multiple_loci_is_not_arbitrarily_picked(self):
        expected = self.root / "expected.tsv"
        expected.write_text("SYNTHETIC1\tderived\n")
        self.write(self.a, [row(), row(pos=2)])
        self.run_cli("--expected-states", str(expected), ok=False)
        self.assert_no_reports()

    def test_callable_row_without_allele_definition_rejected(self):
        self.write(self.a, [row(aa=".")])
        self.run_cli(ok=False)
        self.assert_no_reports()

    def test_unknown_alleles_are_retained_for_noncallable_rows(self):
        self.write(self.a, [row(aa=".", status="unknown")])
        self.write(self.b, [row(aa=".", status="unknown")])
        self.run_cli()
        metrics = self.report("pairwise")[0]
        self.assertEqual(metrics["unverified_allele_keys"], "1")
        self.assertEqual(metrics["callable_both"], "0")

    def test_force_never_overwrites_input_files(self):
        self.a = self.root / "comparison.matrix.tsv"
        self.write(self.a, [row()])
        before = self.a.read_bytes()
        self.run_cli("--force", ok=False)
        self.assertEqual(self.a.read_bytes(), before)

    def test_output_directory_collision_detected_before_writing(self):
        Path(f"{self.out}.summary.tsv").mkdir()
        self.run_cli("--force", ok=False)
        self.assertFalse(Path(f"{self.out}.matrix.tsv").exists())

    def test_old_expected_report_cannot_silently_survive_partial_rerun(self):
        expected = self.root / "expected.tsv"
        expected.write_text("SYNTHETIC1\tderived\n")
        self.run_cli("--expected-states", str(expected))
        self.run_cli("--force", ok=False)


if __name__ == "__main__":
    unittest.main()
