#!/usr/bin/env python3
"""Synthetic resolver regressions; no sample data or external resources.

Run: python3 annotate/test_y_haplo_from_markers.py
Only indexed-FASTA integration checks require pysam; the other checks use stdlib.
"""
from __future__ import annotations

import csv
import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import y_haplo_from_markers as resolver


HERE = Path(__file__).resolve().parent
HAS_PYSAM = importlib.util.find_spec("pysam") is not None
SEQUENCE = "ACGTACGTACGTACGT"
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chrY,length=16>
##ALT=<ID=NON_REF,Description="Unspecified alternate allele">
##ALT=<ID=*,Description="Unspecified alternate allele">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=RefCall,Description="Reference call">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum block depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
"""


def record(pos=1, ref="A", alt="<*>", filt="PASS", end=4, gt="0/0",
           gq="40", dp="10", min_dp="10"):
    info = f"END={end}" if end is not None else "."
    return (f"chrY\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}"
            f"\tGT:GQ:DP:MIN_DP\t{gt}:{gq}:{dp}:{min_dp}\n")


def marker(pos=2, aa="C", der="A", vid="SYNTHETIC1"):
    return resolver.MarkerRow("chrY", pos, vid, aa, der, aa, der,
                              "G-SYNTHETIC", "G2")


class ResolverTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="y-marker-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.vcf = self.root / "input.vcf"
        self.markers = self.root / "markers.vcf"
        self.markers.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chrY,length=16>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chrY\t2\tSYNTHETIC1\tC\tA\t.\t.\tAA=C;HG=G-SYNTHETIC;ISOGG=G2\n"
            "chrY\t3\tSYNTHETIC2\tA\tG\t.\t.\tAA=A;HG=G-SYNTHETIC;ISOGG=G2\n"
        )

    def write_calls(self, *rows):
        self.vcf.write_text(HEADER + "".join(rows))

    def lookup(self, pos, *rows, mode="any"):
        self.write_calls(*rows)
        return resolver.sample_for_pos(pos, *resolver.build_sample_lookup(
            self.vcf, None, "chrY"), filter_mode=mode)

    def classify(self, m, row, **kwargs):
        return resolver.classify_marker(m, row, **kwargs)

    def run_cli(self, *args, ok=True, prefix="result"):
        output = self.root / prefix
        run = subprocess.run(
            [sys.executable, str(HERE / "y_haplo_from_markers.py"),
             "-i", str(self.vcf), "--markers", str(self.markers),
             "-o", str(output), *args], text=True, capture_output=True)
        if ok:
            self.assertEqual(run.returncode, 0, run.stderr)
        else:
            self.assertNotEqual(run.returncode, 0)
        return output, run

    def statuses(self, output):
        with Path(f"{output}.marker_status.tsv").open() as fh:
            return list(csv.reader(fh, delimiter="\t"))

    def make_reference(self):
        path = self.root / "reference.fa"
        path.write_text(f">chrY\n{SEQUENCE}\n")
        Path(f"{path}.fai").write_text("chrY\t16\t6\t16\t17\n")
        return path

    @staticmethod
    def fetch(chrom, start, end):
        if chrom != "chrY":
            raise KeyError(chrom)
        return SEQUENCE[start:end]

    def test_block_not_resolved_from_anchor_without_reference(self):
        row = self.lookup(2, record())
        self.assertEqual(self.classify(marker(), row),
                         ("gvcf_block", "nocall", "reference_required"))

    def test_anchor_can_be_resolved_without_reference(self):
        row = self.lookup(1, record())
        self.assertEqual(self.classify(marker(1, "A", "C"), row),
                         ("gvcf_block", "ancestral", "."))

    def test_block_interior_uses_native_base(self):
        for alt in ("<*>", "<NON_REF>", "."):
            for gt in ("0", "0/0", "0|0"):
                with self.subTest(alt=alt, gt=gt):
                    row = self.lookup(2, record(alt=alt, gt=gt))
                    self.assertTrue(row.is_block)
                    self.assertEqual(self.classify(marker(), row, reference=self.fetch),
                                     ("gvcf_block", "ancestral", "."))
                    self.assertEqual(self.classify(marker(3, "A", "G"), row,
                                                   reference=self.fetch),
                                     ("gvcf_block", "derived", "."))

    def test_block_end_is_inclusive_but_next_position_missing(self):
        self.assertIsNotNone(self.lookup(4, record()))
        self.assertIsNone(self.lookup(5, record()))

    def test_gaps_do_not_inherit_previous_block(self):
        self.assertIsNone(self.lookup(6, record(), record(9, end=12)))

    def test_overlapping_blocks_are_rejected(self):
        with self.assertRaisesRegex(SystemExit, "overlapping reference blocks"):
            self.lookup(3, record(), record(3, "G", end=6))

    def test_backwards_block_interval_is_rejected(self):
        with self.assertRaisesRegex(SystemExit, "END precedes POS"):
            self.lookup(5, record(5, end=4))

    def test_nonreference_genotypes_in_blocks_are_not_base_calls(self):
        for gt, expected in ((".", "nocall"), ("./.", "nocall"),
                             ("0/.", "nocall"), ("0/1", "ambiguous"),
                             ("1", "nocall"), ("1/1", "nocall")):
            with self.subTest(gt=gt):
                row = self.lookup(2, record(gt=gt))
                self.assertEqual(self.classify(marker(), row, reference=self.fetch)[1],
                                 expected)

    def test_concrete_alt_with_nonref_is_not_a_reference_block(self):
        row = self.lookup(2, record(2, "C", "A,<NON_REF>", end=4, gt="1"))
        self.assertFalse(row.is_block)
        self.assertEqual(self.classify(marker(), row)[1], "derived")
        self.assertIsNone(self.lookup(3, record(2, "C", "A,<NON_REF>", end=4, gt="1")))

    def test_symbolic_alt_is_not_a_known_derived_base(self):
        row = self.lookup(2, record(2, "C", "A,<NON_REF>", end=None, gt="2"))
        self.assertEqual(self.classify(marker(), row)[1:],
                         ("nocall", "symbolic_allele"))

    def test_variant_overrides_covering_block(self):
        row = self.lookup(2, record(), record(2, "C", "A", end=None, gt="1"))
        self.assertEqual(self.classify(marker(), row)[1], "derived")

    def test_filtered_variant_does_not_fall_back_to_reference_block(self):
        row = self.lookup(2, record(), record(2, "C", "A", filt="LowQual", end=None, gt="1"),
                          mode="pass")
        self.assertEqual(self.classify(marker(), row, filter_mode="pass"),
                         ("variant", "filtered", "site_filter"))

    def test_passing_duplicate_site_is_preferred_over_filtered_record(self):
        row = self.lookup(2, record(2, "C", "A", filt="LowQual", end=None, gt="1", gq="99"),
                          record(2, "C", "A", end=None, gt="0", gq="30"), mode="pass")
        self.assertEqual(self.classify(marker(), row, filter_mode="pass")[1], "ancestral")

    def test_deepvariant_mode_accepts_unfiltered_blocks_only(self):
        row = self.lookup(2, record(filt="."), mode="deepvariant")
        self.assertEqual(self.classify(marker(), row, reference=self.fetch,
                                      filter_mode="deepvariant")[1], "ancestral")
        row = self.lookup(2, record(2, "C", "A", filt=".", end=None, gt="1"))
        self.assertEqual(self.classify(marker(), row, filter_mode="deepvariant")[1],
                         "filtered")

    def test_quality_thresholds_apply_to_ancestral_and_derived(self):
        for gt in ("0", "1", "0/1"):
            with self.subTest(gt=gt):
                row = self.lookup(2, record(2, "C", "A", end=None, gt=gt, gq="19", dp="2"))
                self.assertEqual(self.classify(marker(), row, min_gq=20, min_dp=3)[1:],
                                 ("filtered", "low_gq;low_dp"))

    def test_threshold_boundaries_and_missing_values(self):
        for gq, dp, status, reason in (("20", "3", "derived", "."),
                                      (".", "3", "filtered", "missing_gq"),
                                      ("20", ".", "filtered", "missing_dp")):
            with self.subTest(gq=gq, dp=dp):
                row = self.lookup(2, record(2, "C", "A", end=None, gt="1", gq=gq, dp=dp))
                self.assertEqual(self.classify(marker(), row, min_gq=20, min_dp=3)[1:],
                                 (status, reason))

    def test_minimum_block_depth_takes_precedence(self):
        row = self.lookup(2, record(dp="30", min_dp="2"))
        self.assertEqual(row.dp, "2")
        self.assertEqual(self.classify(marker(), row, min_dp=3)[1:],
                         ("filtered", "low_dp"))

    def test_block_depth_falls_back_to_dp_only_if_min_dp_absent(self):
        text = record(dp="12").replace(":MIN_DP", "").replace(":40:12:10", ":40:12")
        row = self.lookup(2, text)
        self.assertEqual(row.dp, "12")
        row = self.lookup(2, record(dp="30", min_dp="."))
        self.assertEqual(self.classify(marker(), row, min_dp=3)[1:],
                         ("filtered", "missing_dp"))

    def test_zero_depth_block_cannot_supply_reference_evidence(self):
        row = self.lookup(2, record(dp="0", min_dp="0"))
        self.assertEqual(self.classify(marker(), row, reference=self.fetch)[1:],
                         ("nocall", "zero_depth"))

    def test_mismatched_or_masked_reference_is_not_called(self):
        row = self.lookup(2, record())
        for fetch, reason in ((lambda *_: "T", "reference_mismatch"),
                              (lambda c, s, e: "A" if s == 0 else "N", "noncanonical_base")):
            with self.subTest(reason=reason):
                self.assertEqual(self.classify(marker(), row, reference=fetch)[1:],
                                 ("nocall", reason))

    def test_reference_fetch_failure_is_an_explicit_nocall(self):
        def missing_contig(*args):
            raise KeyError("chrY")
        row = self.lookup(2, record())
        self.assertEqual(self.classify(marker(), row, reference=missing_contig)[1:],
                         ("nocall", "reference_unavailable"))

    def test_interior_indel_marker_not_guessed_from_block(self):
        row = self.lookup(2, record())
        self.assertEqual(self.classify(marker(2, "C", "CA"), row,
                                       reference=self.fetch)[1:],
                         ("nocall", "unsupported_block_marker"))

    def test_variant_only_haploid_and_diploid_calls_remain_supported(self):
        for gt, expected in (("0", "ancestral"), ("0/0", "ancestral"),
                             ("1", "derived"), ("1|1", "derived"),
                             ("0/1", "ambiguous"), ("./.", "nocall")):
            with self.subTest(gt=gt):
                row = self.lookup(2, record(2, "C", "A", end=None, gt=gt))
                self.assertEqual(self.classify(marker(), row)[1], expected)

    def test_reversed_marker_ref_alt_does_not_change_allele_state(self):
        row = self.lookup(2, record(2, "C", "A", end=None, gt="0"))
        m = resolver.MarkerRow("chrY", 2, "SYNTHETIC1", "A", "C", "A", ".", "G2", "G2")
        self.assertEqual(self.classify(m, row)[1], "derived")

    def test_cli_filters_statuses_and_summary_not_just_derived_file(self):
        self.write_calls(record(2, "C", "A", end=None, gt="1", gq="19"),
                         record(3, "G", "A", end=None, gt="0", dp="2"))
        output, _ = self.run_cli("--min-gq", "20", "--min-dp", "3")
        rows = self.statuses(output)
        self.assertEqual([r[12] for r in rows], ["filtered", "filtered"])
        self.assertEqual([r[13] for r in rows], ["low_gq", "low_dp"])
        self.assertTrue(all(len(r)==14 for r in rows))
        self.assertEqual(Path(f"{output}.derived.tsv").read_text(), "")
        summary = Path(f"{output}.summary.txt").read_text()
        for expected in ("  filtered: 2", "  derived: 0", "  resolved_nonmiss: 0"):
            self.assertIn(expected, summary)

    def test_cli_missing_filtered_and_explicit_nocall_are_distinct(self):
        self.write_calls(record(2, "C", "A", end=None, gt="1", filt="LowQual"))
        output, _ = self.run_cli("--site-filter-mode", "pass")
        self.assertEqual([r[11:14] for r in self.statuses(output)],
                         [["variant", "filtered", "site_filter"], ["missing", "nocall", "missing_record"]])
        self.write_calls(record(2, "C", "A", end=None, gt="./."))
        output, _ = self.run_cli(prefix="explicit")
        self.assertEqual(self.statuses(output)[0][11:14], ["variant", "nocall", "missing_gt"])

    def test_gff3_derived_and_unplaceable_output_have_same_schema(self):
        self.markers = self.root / "markers.gff3"
        self.markers.write_text(
            "##gff-version 3\nchrY\tpoint\tsnp\t2\t2\t.\t+\t.\t"
            "ID=SYNTHETIC1;allele_anc=C;allele_der=A;isogg_haplogroup=unknown;yfull_node=G-SYNTHETIC\n")
        self.write_calls(record(2, "C", "A", end=None, gt="1"))
        output, _ = self.run_cli()
        rows = self.statuses(output)
        with Path(f"{output}.derived.tsv").open() as fh:
            self.assertEqual(list(csv.reader(fh, delimiter="\t")), rows)
        with Path(f"{output}.unplaceable_derived.tsv").open() as fh:
            unplaced = list(csv.reader(fh, delimiter="\t"))
        self.assertEqual(unplaced[0], list(resolver.OUTPUT_COLUMNS))
        self.assertEqual(unplaced[1:], rows)

    def test_cli_no_reference_reports_why_interior_blocks_are_unresolved(self):
        self.write_calls(record())
        output, run = self.run_cli()
        self.assertEqual([r[12:14] for r in self.statuses(output)],
                         [["nocall", "reference_required"], ["nocall", "reference_required"]])
        self.assertIn("--reference", run.stderr)

    @unittest.skipUnless(HAS_PYSAM, "indexed FASTA integration requires pysam")
    def test_cli_uses_indexed_reference_at_each_marker(self):
        self.write_calls(record(alt="<NON_REF>", filt="."))
        output, _ = self.run_cli("--reference", str(self.make_reference()),
                                 "--site-filter-mode", "deepvariant")
        self.assertEqual([r[12:14] for r in self.statuses(output)],
                         [["ancestral", "."], ["derived", "."]])

    @unittest.skipUnless(HAS_PYSAM, "BGZF reference integration requires pysam")
    def test_cli_accepts_bgzf_reference_with_its_indexes(self):
        import pysam
        ref = self.make_reference()
        compressed = self.root / "reference.fa.gz"
        pysam.tabix_compress(str(ref), str(compressed))
        pysam.faidx(str(compressed))
        self.write_calls(record())
        output, _ = self.run_cli("--reference", str(compressed))
        self.assertEqual(self.statuses(output)[0][12], "ancestral")

    @unittest.skipUnless(HAS_PYSAM, "VCF regional fetch integration requires pysam")
    def test_indexed_input_matches_sequential_scan(self):
        import pysam
        self.write_calls(record())
        plain, _ = self.run_cli(prefix="plain")
        self.vcf = Path(pysam.tabix_index(str(self.vcf), preset="vcf", keep_original=True))
        indexed, _ = self.run_cli(prefix="indexed")
        self.assertEqual(self.statuses(plain), self.statuses(indexed))

    def test_plain_gzip_reference_is_rejected(self):
        ref = self.root / "reference.fa.gz"
        with gzip.open(ref, "wt") as fh:
            fh.write(f">chrY\n{SEQUENCE}\n")
        Path(f"{ref}.fai").write_text("chrY\t16\t6\t16\t17\n")
        self.write_calls(record())
        _, run = self.run_cli("--reference", str(ref), ok=False)
        self.assertIn("not BGZF", run.stderr)

    def test_missing_fasta_index_is_not_silently_created(self):
        self.write_calls(record())
        ref = self.make_reference()
        Path(f"{ref}.fai").unlink()
        output, run = self.run_cli("--reference", str(ref), ok=False)
        self.assertIn(".fai", run.stderr)
        self.assertFalse(Path(f"{ref}.fai").exists())
        self.assertFalse(Path(f"{output}.marker_status.tsv").exists())

    @unittest.skipUnless(HAS_PYSAM, "malformed FASTA index integration requires pysam")
    def test_malformed_fasta_index_is_not_silently_rebuilt(self):
        self.write_calls(record())
        ref = self.make_reference()
        fai = Path(f"{ref}.fai")
        fai.write_text("not_an_index\n")
        before = fai.read_bytes()
        output, _ = self.run_cli("--reference", str(ref), ok=False)
        self.assertEqual(fai.read_bytes(), before)
        self.assertFalse(Path(f"{output}.marker_status.tsv").exists())

    def test_rejects_negative_threshold(self):
        self.write_calls(record())
        for arg in ("--min-gq", "--min-dp"):
            with self.subTest(arg=arg):
                self.run_cli(arg, "-1", ok=False)

    def test_dual_driver_rejects_bad_threshold_before_work(self):
        for option in ("--min-gq", "--min-dp"):
            run = subprocess.run(["sh", str(HERE / "run_y_dual_liftover_experiment.sh"), option, "-1"],
                                 cwd=self.root, text=True, capture_output=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("nonnegative integers", run.stderr)

    @unittest.skipUnless(HAS_PYSAM and shutil.which("bcftools"),
                         "driver integration requires pysam and bcftools")
    def test_iceman_driver_forwards_reference_and_sample(self):
        import pysam
        self.write_calls(record())
        indexed = Path(pysam.tabix_index(str(self.vcf), preset="vcf", keep_original=True))
        ref = self.make_reference()
        dest = self.root / "output with spaces"
        run = subprocess.run(
            ["sh", str(HERE / "run_iceman_y_compare.sh"),
             "--iceman-gvcf", str(indexed), "--iceman-gvcf-sample", "SYNTHETIC",
             "--markers", str(self.markers), "--reference", str(ref),
             "--out-dir", str(dest), "--skip-vcf", "--skip-solid"],
            cwd=self.root, env={**os.environ, "HOME": str(self.root)},
            text=True, capture_output=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual([r[12] for r in self.statuses(dest / "iceman_gvcf")],
                         ["ancestral", "derived"])

    def test_gzip_content_detected_without_gz_suffix(self):
        self.write_calls(record(2, "C", "A", end=None, gt="1"))
        with gzip.open(self.vcf, "wb") as fh:
            fh.write((HEADER + record(2, "C", "A", end=None, gt="1")).encode())
        output, _ = self.run_cli()
        self.assertEqual(self.statuses(output)[0][12], "derived")

    def test_filtered_rows_do_not_support_ranking_or_create_ancestral_flips(self):
        self.write_calls(record(2, "C", "A", end=None, gt="1", gq="10"),
                         record(3, "G", "A", end=None, gt="0"))
        raw, _ = self.run_cli(prefix="raw")
        filtered, _ = self.run_cli("--min-gq", "20", prefix="filtered")
        ranks = self.root / "ranks.tsv"
        subprocess.run([sys.executable, str(HERE / "y_path_rank.py"),
                        "--input", f"{filtered}.marker_status.tsv", "--out", str(ranks)],
                       text=True, capture_output=True, check=True)
        with ranks.open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        by_candidate = {row["candidate"]: row for row in rows}
        self.assertEqual(by_candidate["G-SYNTHETIC1"]["derived_total"], "0")
        self.assertEqual(by_candidate["G-SYNTHETIC1"]["nocall_total"], "1")
        compared = self.root / "comparison"
        subprocess.run([sys.executable, str(HERE / "compare_marker_status_sets.py"),
                        "--input", f"raw={raw}.marker_status.tsv",
                        "--input", f"filtered={filtered}.marker_status.tsv",
                        "--out-prefix", str(compared)],
                       text=True, capture_output=True, check=True)
        with compared.with_suffix(".matrix.tsv").open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(rows[0]["filtered.status"], "filtered")
        self.assertEqual(compared.with_suffix(".summary.tsv").read_text().splitlines()[-1].split("\t")[-1], "0")

    def test_filtered_defining_marker_does_not_deepen_tree_placement(self):
        self.write_calls(record(2, "C", "A", end=None, gt="1", gq="10"),
                         record(3, "G", "A", end=None, gt="0"))
        output, _ = self.run_cli("--min-gq", "20")
        tree = self.root / "tree.tsv"
        tree.write_text("node\tparent\tstatus\tdefining\tequivalent\tevidence\n"
                        "G-ROOT\t-\tpublished\tSYNTHETIC2\t-\tsynthetic\n"
                        "G-CHILD\tG-ROOT\tpublished\tSYNTHETIC1\t-\tsynthetic\n")
        tree_out = self.root / "tree_calls.tsv"
        subprocess.run([sys.executable, str(HERE / "y_path_rank.py"),
                        "--input", f"{output}.marker_status.tsv", "--out", str(self.root / "ranks.tsv"),
                        "--tree", str(tree), "--tree-out", str(tree_out)],
                       capture_output=True, text=True, check=True)
        with tree_out.open() as fh:
            nodes = {r["node"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(nodes["G-ROOT"]["terminal"], "yes")
        self.assertEqual(nodes["G-CHILD"]["terminal"], "no")
        self.assertEqual(nodes["G-CHILD"]["def_nocall"], "1")

    def test_legacy_13_column_tables_remain_readable(self):
        self.write_calls(record(2, "C", "A", end=None, gt="1"))
        output, _ = self.run_cli()
        legacy = self.root / "legacy.tsv"
        legacy.write_text("".join("\t".join(r[:13]) + "\n" for r in self.statuses(output)))
        ranking = self.root / "legacy_ranks.tsv"
        subprocess.run([sys.executable, str(HERE / "y_path_rank.py"),
                        "--input", str(legacy), "--out", str(ranking)],
                       capture_output=True, text=True, check=True)
        with ranking.open() as fh:
            candidates = {r["candidate"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(candidates["G-SYNTHETIC1"]["derived_total"], "1")
        comparison = self.root / "legacy_comparison"
        subprocess.run([sys.executable, str(HERE / "compare_marker_status_sets.py"),
                        "--input", f"new={output}.marker_status.tsv", "--input", f"old={legacy}",
                        "--out-prefix", str(comparison)], capture_output=True, text=True, check=True)
        self.assertEqual(comparison.with_suffix(".summary.tsv").read_text().splitlines()[-1].split("\t")[2:],
                         ["2", "0", "0"])


if __name__ == "__main__":
    unittest.main()
