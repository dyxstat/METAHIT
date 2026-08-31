#!/usr/bin/env python3
"""Run the packaged Binning_refiner and preserve METAHICT's output contract.

This file contains only integration code.  The refinement algorithm is supplied
by the upstream Binning_refiner 1.4.3 package pinned in
installation/pip-requirements.txt.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-size", required=True, type=int, help="Minimum bin size in bp")
    parser.add_argument("bin_sets", nargs="+", type=Path)
    return parser.parse_args()


def minimum_kbp_for_upstream(minimum_bp: int) -> int:
    """Translate METAHICT's base-pair cutoff to upstream's integer KiB unit."""
    if minimum_bp < 0:
        raise ValueError("minimum bin size cannot be negative")
    return math.ceil(minimum_bp / 1024)


def main() -> int:
    args = parse_args()
    if len(args.bin_sets) < 2:
        raise SystemExit("Binning_refiner requires at least two bin sets")

    executable = os.environ.get("METAHICT_BINNING_REFINER_BIN") or shutil.which("Binning_refiner")
    if not executable:
        raise SystemExit(
            "Binning_refiner was not found. Run ./metahict install to install "
            "the pinned upstream Binning_refiner 1.4.3 package."
        )

    bin_sets = [path.resolve(strict=True) for path in args.bin_sets]
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        shutil.rmtree(output)

    with tempfile.TemporaryDirectory(prefix="metahict-binning-refiner-") as tmp:
        workspace = Path(tmp)
        inputs = workspace / "inputs"
        inputs.mkdir()
        for index, bin_set in enumerate(bin_sets, start=1):
            (inputs / f"set{index}").symlink_to(bin_set, target_is_directory=True)

        # Upstream accepts a zero-KiB cutoff. Preserve METAHICT's documented
        # zero-bp test setting instead of silently raising it to 1 KiB.
        minimum_kbp = minimum_kbp_for_upstream(args.min_size)
        subprocess.run(
            [str(executable), "-i", str(inputs), "-p", "METAHICT", "-m", str(minimum_kbp)],
            cwd=workspace,
            check=True,
        )
        upstream_output = workspace / "METAHICT_Binning_refiner_outputs"
        refined_bins = upstream_output / "METAHICT_refined_bins"
        if not refined_bins.is_dir():
            raise SystemExit(f"Binning_refiner did not create its expected output: {refined_bins}")

        output.mkdir()
        shutil.copytree(refined_bins, output / "Refined")
        for report in upstream_output.glob("METAHICT_*.txt"):
            shutil.copy2(report, output / report.name)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
