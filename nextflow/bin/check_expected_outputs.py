#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check expected output files for a METAHICT Nextflow run."
    )
    parser.add_argument("--root", required=True, help="Repository or output root")
    parser.add_argument("--manifest", required=True, help="TSV expected-output manifest")
    return parser.parse_args()


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def main():
    args = parse_args()
    root = Path(args.root).resolve()
    manifest = Path(args.manifest).resolve()
    failures = []

    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rel_path = row["path"]
            uses_glob = any(character in rel_path for character in "*?[")
            paths = sorted(root.glob(rel_path)) if uses_glob else [root / rel_path]
            required = truthy(row.get("required", "1"))
            non_empty = truthy(row.get("non_empty", "0"))
            description = row.get("description", "")

            existing = [path for path in paths if path.exists()]
            if required and not existing:
                failures.append(f"missing\t{rel_path}\t{description}")
                continue
            if non_empty:
                for path in existing:
                    if path.stat().st_size == 0:
                        failures.append(
                            f"empty\t{path.relative_to(root)}\t{description}"
                        )

    if failures:
        print("Expected-output check failed:")
        for failure in failures:
            print(failure)
        raise SystemExit(1)

    print(f"Expected-output check passed: {manifest}")


if __name__ == "__main__":
    main()
