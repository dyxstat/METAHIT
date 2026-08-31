#!/usr/bin/env python3
"""Fail CI if compiled tools or known copied implementations enter modules/."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
MODULES = ROOT / "modules"
FORBIDDEN_NAMES = {
    "jgi_summarize_bam_contig_depths",
    "refine_bin_intersections.py",
}
EXECUTABLE_MAGICS = (
    b"\x7fELF",
    b"MZ",
    b"\xcf\xfa\xed\xfe",
    b"\xce\xfa\xed\xfe",
    b"\xfe\xed\xfa\xcf",
    b"\xfe\xed\xfa\xce",
)


def is_macos_metadata_path(path: Path) -> bool:
    return any(
        part == "__MACOSX" or part == ".DS_Store" or part.startswith("._")
        for part in path.parts
    )


def main() -> int:
    violations: list[str] = []
    for path in MODULES.rglob("*"):
        if not path.is_file() or is_macos_metadata_path(path):
            continue
        if path.name in FORBIDDEN_NAMES:
            violations.append(f"known third-party tool copy: {path.relative_to(ROOT)}")
        with path.open("rb") as handle:
            magic = handle.read(4)
        if any(magic.startswith(signature) for signature in EXECUTABLE_MAGICS):
            violations.append(f"compiled executable in modules/: {path.relative_to(ROOT)}")
        if path.stat().st_size > 1024 * 1024:
            violations.append(f"file larger than 1 MiB in modules/: {path.relative_to(ROOT)}")

    if violations:
        print("Vendored-tool policy violations:", file=sys.stderr)
        for violation in violations:
            print(f"- {violation}", file=sys.stderr)
        return 1
    print("No vendored compiled tools or known copied implementations found in modules/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
