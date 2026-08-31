#!/usr/bin/env python3
"""Plot ranked CheckM2 and N50 statistics before and after reassembly."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


class PlotError(RuntimeError):
    """Expected statistics or plotting failure."""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot METAHICT reassembly quality rankings"
    )
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("min_completeness", type=float)
    parser.add_argument("max_contamination", type=float)
    parser.add_argument("stats", nargs="+", type=Path)
    return parser


def validate_thresholds(min_completeness: float, max_contamination: float) -> None:
    if not 0 <= min_completeness <= 100:
        raise PlotError("Minimum completeness must be between 0 and 100")
    if not 0 <= max_contamination <= 100:
        raise PlotError("Maximum contamination must be between 0 and 100")


def load_scores(
    path: Path, min_completeness: float, max_contamination: float
) -> tuple[list[int], list[float], list[float]]:
    if not path.is_file():
        raise PlotError(f"Statistics file not found: {path}")
    n50_values: list[int] = []
    completions: list[float] = []
    contaminations: list[float] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) >= 2 and fields[0].lower() == "bin" and "compl" in fields[1].lower():
                continue
            if len(fields) < 6:
                raise PlotError(f"Malformed statistics row {path}:{line_number}")
            try:
                completeness = float(fields[1])
                contamination = float(fields[2])
                n50 = int(float(fields[5]))
            except ValueError as error:
                raise PlotError(
                    f"Non-numeric CheckM2 value at {path}:{line_number}"
                ) from error
            if completeness >= min_completeness and contamination <= max_contamination:
                n50_values.append(n50)
                completions.append(completeness)
                contaminations.append(contamination)
    n50_values.sort(reverse=True)
    completions.sort(reverse=True)
    contaminations.sort()
    return n50_values, completions, contaminations


def plot_rankings(args: argparse.Namespace) -> None:
    validate_thresholds(args.min_completeness, args.max_contamination)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    datasets = [
        (path.stem, *load_scores(path, args.min_completeness, args.max_contamination))
        for path in args.stats
    ]

    plt.style.use("ggplot")
    figure, (n50_axis, completion_axis, contamination_axis) = plt.subplots(
        1, 3, figsize=(18, 8)
    )
    colors = plt.get_cmap("tab20").colors
    maximum_rank = 1
    for index, (label, n50_values, completions, contaminations) in enumerate(datasets):
        color = colors[index % len(colors)]
        for axis, values in (
            (n50_axis, n50_values),
            (completion_axis, completions),
            (contamination_axis, contaminations),
        ):
            axis.plot(
                range(1, len(values) + 1),
                values,
                linewidth=2.5,
                color=color,
                label=label,
            )
            maximum_rank = max(maximum_rank, len(values))

    for axis in (n50_axis, completion_axis, contamination_axis):
        axis.set_facecolor("white")
        axis.grid(True, linestyle="--", linewidth=0.5, color="black", alpha=0.3)
        axis.set_xlim(0.5, maximum_rank + 0.5)
        axis.tick_params(labelsize=11)

    n50_axis.set_ylim(bottom=0)
    n50_axis.set_title("Bin N50 ranking", fontsize=18)
    n50_axis.set_xlabel("Descending N50 rank", fontsize=13)
    n50_axis.set_ylabel("Bin N50 (bp)", fontsize=13)

    completion_axis.set_ylim(args.min_completeness, 105)
    completion_axis.yaxis.set_major_formatter(PercentFormatter(xmax=100))
    completion_axis.set_title("Bin completeness ranking", fontsize=18)
    completion_axis.set_xlabel("Descending completeness rank", fontsize=13)
    completion_axis.set_ylabel("Estimated bin completeness", fontsize=13)

    contamination_margin = max(0.5, args.max_contamination * 0.05)
    contamination_axis.set_ylim(0, args.max_contamination + contamination_margin)
    contamination_axis.yaxis.set_major_formatter(PercentFormatter(xmax=100))
    contamination_axis.set_title("Bin contamination ranking", fontsize=18)
    contamination_axis.set_xlabel("Ascending contamination rank", fontsize=13)
    contamination_axis.set_ylabel("Estimated bin contamination", fontsize=13)
    contamination_axis.legend(loc="best", fontsize=9)

    figure.tight_layout(w_pad=3)
    png = args.output_directory / "reassembly_results.png"
    eps = args.output_directory / "reassembly_results.eps"
    figure.savefig(png, format="png", dpi=300)
    figure.savefig(eps, format="eps", dpi=300)
    plt.close(figure)
    print(f"[PASS] Saved {png} and {eps}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        plot_rankings(build_parser().parse_args(argv))
    except (OSError, PlotError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
