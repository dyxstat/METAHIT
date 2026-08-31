#!/usr/bin/env python3
"""Unit checks for METAHICT's upstream Binning_refiner adapter."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ADAPTER = (
    PROJECT_ROOT
    / "modules"
    / "6_binning"
    / "run_binning_refiner.py"
)
SPEC = importlib.util.spec_from_file_location("metahict_binning_refiner_adapter", ADAPTER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load adapter: {ADAPTER}")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MinimumSizeConversionTest(unittest.TestCase):
    def test_zero_is_preserved(self) -> None:
        self.assertEqual(MODULE.minimum_kbp_for_upstream(0), 0)

    def test_exact_kib_is_preserved(self) -> None:
        self.assertEqual(MODULE.minimum_kbp_for_upstream(524_288), 512)

    def test_partial_kib_rounds_up(self) -> None:
        self.assertEqual(MODULE.minimum_kbp_for_upstream(1_025), 2)

    def test_negative_size_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.minimum_kbp_for_upstream(-1)


if __name__ == "__main__":
    unittest.main()
