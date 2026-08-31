#!/usr/bin/env python3
"""Compare accepted and candidate METAHICT scientific result trees."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys


LEGACY_RESULT_PREFIXES = (
    ("10_MGE/intermediates/alignment/", ("10_MGE/alignment/", "10_MGE/mge_alignment/mge_alignment/")),
    ("10_MGE/intermediates/contact/", ("10_MGE/contact/", "10_MGE/mge_contact/mge_contact/")),
    ("10_MGE/", ("10_MGE/results/", "10_MGE/mge/")),
    ("2_assembly/", ("2_assembly/assembly/",)),
    ("3_alignment/", ("3_alignment/alignment/",)),
    ("4_coverage/", ("4_coverage/coverage/",)),
    ("5_contact/", ("5_contact/contact/",)),
    ("6_binning/", ("6_binning/binning/",)),
    ("7_reassembly/", ("7_reassembly/reassembly/",)),
    ("8_scaffolding/", ("8_scaffolding/scaffolding/",)),
    ("9_annotation/", ("9_annotation/annotation/",)),
)


@dataclass
class Result:
    path: str
    kind: str
    status: str
    detail: str


def compatible_result_path(root: Path, relative: str) -> Path:
    """Resolve the current result contract, falling back to the legacy layout."""
    current = root / relative
    if current.exists():
        return current
    for current_prefix, legacy_prefixes in LEGACY_RESULT_PREFIXES:
        if relative.startswith(current_prefix):
            for legacy_prefix in legacy_prefixes:
                legacy = root / (legacy_prefix + relative[len(current_prefix) :])
                if legacy.exists():
                    return legacy
    return current


def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    sequence: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence).upper()
                identifier = line[1:].split()[0]
                if identifier in records:
                    raise ValueError(f"Duplicate FASTA identifier in {path}: {identifier}")
                sequence = []
            elif identifier is None:
                raise ValueError(f"Sequence occurs before a FASTA header in {path}")
            else:
                sequence.append(line)
    if identifier is not None:
        records[identifier] = "".join(sequence).upper()
    return records


def compare_fasta(baseline: Path, candidate: Path) -> tuple[bool, str]:
    left, right = fasta_records(baseline), fasta_records(candidate)
    missing = sorted(left.keys() - right.keys())
    extra = sorted(right.keys() - left.keys())
    changed = sorted(name for name in left.keys() & right.keys() if left[name] != right[name])
    passed = not missing and not extra and not changed
    detail = f"records={len(left)}; missing={len(missing)}; extra={len(extra)}; sequence_changes={len(changed)}"
    if not passed:
        detail += f"; examples={(missing + extra + changed)[:5]}"
    return passed, detail


def table_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first else ","
        rows = list(csv.reader(handle, delimiter=delimiter))
    if not rows:
        return [], []
    return rows[0], sorted(rows[1:])


def numeric(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def compare_table(baseline: Path, candidate: Path, absolute: float, relative: float) -> tuple[bool, str]:
    left_header, left_rows = table_rows(baseline)
    right_header, right_rows = table_rows(candidate)
    if left_header != right_header:
        return False, "column headers differ"
    if len(left_rows) != len(right_rows):
        return False, f"row counts differ: baseline={len(left_rows)}, candidate={len(right_rows)}"
    differences = 0
    maximum_delta = 0.0
    for left, right in zip(left_rows, right_rows):
        if len(left) != len(right):
            differences += 1
            continue
        for old, new in zip(left, right):
            old_number, new_number = numeric(old), numeric(new)
            if old_number is None or new_number is None:
                differences += old != new
                continue
            delta = abs(old_number - new_number)
            maximum_delta = max(maximum_delta, delta)
            scale = max(abs(old_number), abs(new_number))
            differences += delta > absolute + relative * scale
    return differences == 0, f"rows={len(left_rows)}; differing_cells={differences}; max_numeric_delta={maximum_delta:g}"


def fasta_files(directory: Path) -> dict[str, Path]:
    return {
        str(path.relative_to(directory)): path
        for path in directory.rglob("*")
        if path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna"}
    }


def compare_fasta_directory(baseline: Path, candidate: Path) -> tuple[bool, str]:
    left, right = fasta_files(baseline), fasta_files(candidate)
    if left.keys() != right.keys():
        missing = sorted(left.keys() - right.keys())
        extra = sorted(right.keys() - left.keys())
        return False, f"FASTA file sets differ; missing={missing[:5]}; extra={extra[:5]}"
    changed = []
    for relative in sorted(left):
        passed, _ = compare_fasta(left[relative], right[relative])
        if not passed:
            changed.append(relative)
    return not changed, f"files={len(left)}; changed={len(changed)}; examples={changed[:5]}"


def compare_sparse(baseline: Path, candidate: Path, absolute: float, relative: float) -> tuple[bool, str]:
    try:
        import numpy as np
        from scipy.sparse import load_npz
    except ImportError as error:
        raise RuntimeError("numpy and scipy are required to compare sparse matrices") from error
    left, right = load_npz(baseline).tocsr(), load_npz(candidate).tocsr()
    if left.shape != right.shape:
        return False, f"matrix shapes differ: baseline={left.shape}, candidate={right.shape}"
    difference = (left - right).tocoo()
    if difference.nnz == 0:
        return True, f"shape={left.shape}; differing_entries=0"
    deltas = np.abs(difference.data)
    scale = max(float(np.max(np.abs(left.data))) if left.nnz else 0.0,
                float(np.max(np.abs(right.data))) if right.nnz else 0.0)
    threshold = absolute + relative * scale
    count = int(np.count_nonzero(deltas > threshold))
    return count == 0, f"shape={left.shape}; entries_over_tolerance={count}; max_delta={float(deltas.max()):g}"


def compare_entry(kind: str, baseline: Path, candidate: Path, absolute: float, relative: float) -> tuple[bool, str]:
    if kind == "fasta":
        return compare_fasta(baseline, candidate)
    if kind == "fasta_dir":
        return compare_fasta_directory(baseline, candidate)
    if kind == "table":
        return compare_table(baseline, candidate, absolute, relative)
    if kind == "sparse_npz":
        return compare_sparse(baseline, candidate, absolute, relative)
    if kind == "exact":
        passed = baseline.read_bytes() == candidate.read_bytes()
        return passed, "byte-identical" if passed else "file bytes differ"
    raise ValueError(f"Unsupported comparison kind: {kind}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    results: list[Result] = []
    with args.manifest.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            relative, kind = row["path"], row["kind"]
            baseline = compatible_result_path(args.baseline, relative)
            candidate = compatible_result_path(args.candidate, relative)
            if not baseline.exists() or not candidate.exists():
                missing = [label for label, path in (("baseline", baseline), ("candidate", candidate)) if not path.exists()]
                results.append(Result(relative, kind, "FAIL", "missing from " + ", ".join(missing)))
                continue
            try:
                passed, detail = compare_entry(
                    kind,
                    baseline,
                    candidate,
                    float(row.get("absolute_tolerance") or 0),
                    float(row.get("relative_tolerance") or 0),
                )
            except (OSError, RuntimeError, ValueError) as error:
                passed, detail = False, str(error)
            results.append(Result(relative, kind, "PASS" if passed else "FAIL", detail))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "kind", "status", "detail"))
        for result in results:
            writer.writerow((result.path, result.kind, result.status, result.detail))
    failures = [result for result in results if result.status == "FAIL"]
    for result in results:
        print(f"[{result.status}] {result.path}: {result.detail}")
    print(f"Scientific comparison: {len(results) - len(failures)} passed, {len(failures)} failed; report={args.report}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
