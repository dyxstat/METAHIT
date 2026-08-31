#!/usr/bin/env python3
"""Create tiny, valid read files and a samplesheet for Nextflow stub tests."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


READ_NAMES = (
    "sg_R1.fastq.gz",
    "sg_R2.fastq.gz",
    "hic_R1.fastq.gz",
    "hic_R2.fastq.gz",
)
LONG_READ_NAME = "long_reads.fastq.gz"


def create_stub_inputs(output_directory: Path) -> Path:
    output_directory = output_directory.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    for index, name in enumerate(READ_NAMES, start=1):
        with gzip.open(output_directory / name, "wt") as handle:
            handle.write(f"@stub_read_{index}\nACGTACGT\n+\nIIIIIIII\n")
    with gzip.open(output_directory / LONG_READ_NAME, "wt") as handle:
        handle.write("@stub_long_read\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")

    samplesheet = output_directory / "samplesheet.csv"
    with samplesheet.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ("sample", "sg1", "sg2", "hic1", "hic2", "enzyme", "long_read_type")
        )
        writer.writerow(
            (
                "example_dataset",
                *(str((output_directory / name).resolve()) for name in READ_NAMES),
                "Sau3AI,MluCI",
                "",
            )
        )
        writer.writerow(
            (
                "long_read_example",
                str((output_directory / LONG_READ_NAME).resolve()),
                "",
                str((output_directory / "hic_R1.fastq.gz").resolve()),
                str((output_directory / "hic_R2.fastq.gz").resolve()),
                "Sau3AI,MluCI",
                "nano-hq",
            )
        )
    return samplesheet


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory in which to create the temporary reads and samplesheet",
    )
    args = parser.parse_args()
    print(create_stub_inputs(args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
