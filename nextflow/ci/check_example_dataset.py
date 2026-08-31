#!/usr/bin/env python3
"""Validate the bundled paired FASTQ example without scientific dependencies."""

from __future__ import annotations

from contextlib import ExitStack
import csv
import gzip
from itertools import zip_longest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
DATASET = ROOT / "example_dataset"
SAMPLESHEET = ROOT / "nextflow" / "assets" / "example_dataset_samplesheet.csv"
EXPECTED_PAIRS = {"sg": 199_981, "hic": 92_384}


def canonical_name(header: str) -> str:
    name = header[1:].split(maxsplit=1)[0]
    return name[:-2] if name.endswith(("/1", "/2")) else name


def fastq_records(handle, path: Path):
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        if not sequence or not separator or not quality:
            raise ValueError(f"Truncated FASTQ record in {path}")
        if not header.startswith("@") or not separator.startswith("+"):
            raise ValueError(f"Invalid FASTQ record in {path}: {header.rstrip()}")
        if len(sequence.rstrip()) != len(quality.rstrip()):
            raise ValueError(f"Sequence/quality length mismatch in {path}: {header.rstrip()}")
        yield header


def validate_pair(library: str) -> int:
    forward = DATASET / f"{library}_R1.fastq.gz"
    reverse = DATASET / f"{library}_R2.fastq.gz"
    with ExitStack() as stack:
        r1 = stack.enter_context(gzip.open(forward, "rt"))
        r2 = stack.enter_context(gzip.open(reverse, "rt"))
        count = 0
        for count, records in enumerate(
            zip_longest(fastq_records(r1, forward), fastq_records(r2, reverse)),
            start=1,
        ):
            header1, header2 = records
            if header1 is None or header2 is None:
                raise ValueError(f"Unequal {library} R1/R2 record counts")
            if canonical_name(header1) != canonical_name(header2):
                raise ValueError(
                    f"Unsynchronized {library} pair {count}: "
                    f"{header1.rstrip()} != {header2.rstrip()}"
                )
    if count != EXPECTED_PAIRS[library]:
        raise ValueError(
            f"Unexpected {library} pair count: {count}; expected {EXPECTED_PAIRS[library]}"
        )
    return count


def validate_samplesheet() -> None:
    with SAMPLESHEET.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1 or rows[0].get("sample") != "example_dataset":
        raise ValueError("Example samplesheet must contain one 'example_dataset' row")
    expected_files = {
        "sg1": "sg_R1.fastq.gz",
        "sg2": "sg_R2.fastq.gz",
        "hic1": "hic_R1.fastq.gz",
        "hic2": "hic_R2.fastq.gz",
    }
    for field, filename in expected_files.items():
        configured = (ROOT / rows[0][field]).resolve()
        expected = (DATASET / filename).resolve()
        if configured != expected or not configured.is_file():
            raise ValueError(f"Invalid example samplesheet path for {field}: {configured}")
    if rows[0].get("enzyme") != "Sau3AI,MluCI":
        raise ValueError("Example samplesheet enzyme must be Sau3AI,MluCI")


def main() -> int:
    try:
        validate_samplesheet()
        shotgun = validate_pair("sg")
        hic = validate_pair("hic")
    except (OSError, UnicodeError, ValueError) as error:
        print(f"[ERROR] Example dataset validation failed: {error}", file=sys.stderr)
        return 1
    print(
        "[PASS] Example dataset: "
        f"{shotgun:,} paired shotgun reads and {hic:,} paired Hi-C reads"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
