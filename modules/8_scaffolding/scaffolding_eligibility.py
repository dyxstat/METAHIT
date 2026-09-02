#!/usr/bin/env python3
"""Dependency-free FASTA eligibility checks for the scaffolding module."""

from __future__ import annotations

import csv
from pathlib import Path


STATUS_FIELDS = (
    "bin",
    "status",
    "reason",
    "total_contigs",
    "eligible_contigs",
    "longest_contig_bp",
    "min_contig_length",
)


def fasta_sequence_lengths(fasta: str | Path) -> dict[str, int]:
    """Return sequence lengths keyed by the first token of each FASTA header."""
    lengths: dict[str, int] = {}
    name: str | None = None
    length = 0
    with Path(fasta).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if name is not None:
                    if name in lengths:
                        raise ValueError(f"Duplicate FASTA sequence identifier: {name}")
                    lengths[name] = length
                fields = line[1:].strip().split()
                if not fields:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                name = fields[0]
                length = 0
            else:
                sequence = line.strip()
                if sequence and name is None:
                    raise ValueError(
                        f"FASTA sequence occurs before the first header at line {line_number}"
                    )
                length += len(sequence)
    if name is not None:
        if name in lengths:
            raise ValueError(f"Duplicate FASTA sequence identifier: {name}")
        lengths[name] = length
    if not lengths:
        raise ValueError("Input FASTA contains no sequences")
    return lengths


def assess_scaffolding_eligibility(
    fasta: str | Path, min_contig_length: int
) -> dict[str, int | str]:
    """Summarize whether a bin has at least two contigs eligible for joining."""
    lengths = fasta_sequence_lengths(fasta)
    eligible = sum(length >= min_contig_length for length in lengths.values())
    reason = ""
    if eligible == 0:
        reason = "no_contigs_meet_minimum_length"
    elif eligible == 1:
        reason = "fewer_than_two_eligible_contigs"
    return {
        "total_contigs": len(lengths),
        "eligible_contigs": eligible,
        "longest_contig_bp": max(lengths.values()),
        "min_contig_length": min_contig_length,
        "reason": reason,
    }


def write_scaffolding_status(
    output: str | Path,
    input_fasta: str | Path,
    status: str,
    assessment: dict[str, int | str],
) -> None:
    """Write one machine-readable status row for a scaffolding attempt."""
    row = {
        "bin": Path(input_fasta).name,
        "status": status,
        "reason": assessment.get("reason", ""),
        "total_contigs": assessment["total_contigs"],
        "eligible_contigs": assessment["eligible_contigs"],
        "longest_contig_bp": assessment["longest_contig_bp"],
        "min_contig_length": assessment["min_contig_length"],
    }
    with Path(output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
