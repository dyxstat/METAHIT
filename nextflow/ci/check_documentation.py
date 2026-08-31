#!/usr/bin/env python3
"""Check local Markdown links and fenced code blocks in user documentation."""

from __future__ import annotations

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[2]


def is_macos_metadata_path(path: Path) -> bool:
    return any(
        part == "__MACOSX" or part == ".DS_Store" or part.startswith("._")
        for part in path.parts
    )


DOCUMENTS = tuple(
    sorted(
        {
            ROOT / "README.md",
            *[
                path
                for path in (ROOT / "docs").rglob("*.md")
                if not is_macos_metadata_path(path)
            ],
            *[
                path
                for path in (ROOT / "nextflow").glob("*.md")
                if not is_macos_metadata_path(path)
            ],
            ROOT / "packaging" / "bioconda" / "README.md",
        }
    )
)


def main() -> int:
    errors: list[str] = []
    for path in DOCUMENTS:
        text = path.read_text()
        if text.count("```") % 2:
            errors.append(f"Unbalanced code fence: {path}")
        for target in re.findall(r"\[[^]]+\]\(([^)]+)\)", text):
            if "://" in target or target.startswith("#"):
                continue
            local_target = target.split("#", 1)[0]
            if local_target and not (path.parent / local_target).resolve().exists():
                errors.append(f"Broken local link in {path}: {target}")
    if errors:
        for error in errors:
            print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    print("Documentation links and code fences passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
