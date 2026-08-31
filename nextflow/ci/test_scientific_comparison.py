#!/usr/bin/env python3
"""Tests for the scientific result comparison utility."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "nextflow" / "bin" / "compare_scientific_outputs.py"
SPEC = importlib.util.spec_from_file_location("scientific_comparison", SCRIPT)
COMPARISON = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = COMPARISON
SPEC.loader.exec_module(COMPARISON)


class ScientificComparisonTest(unittest.TestCase):
    def test_current_manifest_path_resolves_old_nested_result(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            legacy = root / "2_assembly" / "assembly" / "final_assembly.fasta"
            legacy.parent.mkdir(parents=True)
            legacy.write_text(">contig\nAAAA\n")

            observed = COMPARISON.compatible_result_path(
                root, "2_assembly/final_assembly.fasta"
            )

            self.assertEqual(observed, legacy)

    def test_current_path_takes_precedence_over_legacy_result(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            current = root / "10_MGE" / "candidate_mge_host_pairs_zscore_filtered.tsv"
            legacy = root / "10_MGE" / "results" / "candidate_mge_host_pairs_zscore_filtered.tsv"
            current.parent.mkdir(parents=True)
            legacy.parent.mkdir(parents=True)
            current.touch()
            legacy.touch()

            observed = COMPARISON.compatible_result_path(
                root, "10_MGE/candidate_mge_host_pairs_zscore_filtered.tsv"
            )

            self.assertEqual(observed, current)

    def test_flat_mge_path_resolves_previous_results_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            previous = root / "10_MGE" / "results" / "sequence_topology.tsv"
            previous.parent.mkdir(parents=True)
            previous.touch()

            observed = COMPARISON.compatible_result_path(
                root, "10_MGE/sequence_topology.tsv"
            )

            self.assertEqual(observed, previous)

    def test_fasta_comparison_ignores_record_and_line_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old, new = root / "old.fa", root / "new.fa"
            old.write_text(">a description\nAAAA\n>b\nCC\nCC\n")
            new.write_text(">b\nCCCC\n>a other description\nAA\nAA\n")
            passed, detail = COMPARISON.compare_fasta(old, new)
            self.assertTrue(passed, detail)

    def test_table_comparison_applies_numeric_tolerance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old, new = root / "old.tsv", root / "new.tsv"
            old.write_text("name\tvalue\na\t1.000000\n")
            new.write_text("name\tvalue\na\t1.000001\n")
            self.assertTrue(COMPARISON.compare_table(old, new, 1e-5, 0)[0])
            self.assertFalse(COMPARISON.compare_table(old, new, 1e-8, 0)[0])


if __name__ == "__main__":
    unittest.main()
