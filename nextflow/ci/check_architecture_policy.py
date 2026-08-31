#!/usr/bin/env python3
"""Enforce METAHICT's Nextflow/Python orchestration boundary."""

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


PRIMARY_WORKFLOW = (
    ROOT / "nextflow" / "main_dsl2.nf",
    ROOT / "nextflow" / "modules" / "local" / "metahict_modules.nf",
)
BEGINNER_DOCS = (
    ROOT / "README.md",
    ROOT / "docs" / "quickstart.md",
    ROOT / "docs" / "installation.md",
    ROOT / "docs" / "test_dataset.md",
)
def main() -> int:
    errors: list[str] = []
    for path in PRIMARY_WORKFLOW:
        text = path.read_text()
        if re.search(r"\bbash\b", text) or re.search(r"\b[\w./-]+\.sh\b", text):
            errors.append(f"Primary workflow must call named Python interfaces, not shell scripts: {path}")
        if re.search(r"\b(?:for|while)\s*\(", text):
            errors.append(
                f"Nextflow 26 strict syntax requires collection operators instead of imperative loops: {path}"
            )
        method_names = re.findall(r"^def\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(", text, re.MULTILINE)
        duplicate_methods = sorted(
            name for name in set(method_names) if method_names.count(name) > 1
        )
        if duplicate_methods:
            errors.append(
                "Nextflow 26 strict syntax does not support these repeated dynamic "
                f"method names in {path}: {', '.join(duplicate_methods)}"
            )

    legacy_command = re.compile(r"^\s*bash\s+\S+\.sh(?:\s|$)", re.MULTILINE)
    for path in BEGINNER_DOCS:
        if legacy_command.search(path.read_text()):
            errors.append(f"Beginner documentation reintroduced a 'bash script.sh' command: {path}")

    unexpected_shell = sorted(
        path
        for source_root in (ROOT / "modules", ROOT / "installation")
        for path in source_root.rglob("*.sh")
        if not is_macos_metadata_path(path)
    )
    if unexpected_shell:
        errors.append(
            "Scientific modules and installation logic must not use shell scripts: "
            + ", ".join(str(path.relative_to(ROOT)) for path in unexpected_shell)
        )

    unsafe_python = re.compile(r"\bos\.system\s*\(|\bshell\s*=\s*True\b")
    for path in (
        candidate
        for candidate in (ROOT / "modules").rglob("*.py")
        if not is_macos_metadata_path(candidate)
    ):
        if unsafe_python.search(path.read_text()):
            errors.append(f"Scientific Python must use checked argument-vector subprocesses: {path}")

    for module in (ROOT / "modules").iterdir():
        if not module.is_dir() or module.name == "__pycache__":
            continue
        nested = [
            path
            for path in module.rglob("*")
            if path.is_dir()
            and path.name != "__pycache__"
            and not is_macos_metadata_path(path)
        ]
        if nested:
            errors.append(
                f"Numbered module must be flat; found subdirectories: "
                f"{', '.join(str(path.relative_to(ROOT)) for path in nested)}"
            )

    workflow_text = (ROOT / "nextflow" / "modules" / "local" / "metahict_modules.nf").read_text()
    for stage in (
        "modules/1_preprocessing/preprocessing.py",
        "modules/2_assembly/assembly.py",
        "modules/3_alignment/alignment.py",
        "modules/4_coverage/coverage.py",
        "modules/5_contact/contact.py",
        "modules/8_scaffolding/scaffolding.py",
        "modules/9_annotation/annotation.py",
        "run_binning.py",
        "run_reassembly.py",
        "modules/10_MGE/mge.py",
    ):
        if stage not in workflow_text:
            errors.append(f"Nextflow does not call the converted Python stage directly: {stage}")

    analytical_blocks = re.findall(
        r"^process (?!CONDA_BUNDLE)\w+ \{.*?(?=^process |\Z)",
        workflow_text,
        re.MULTILINE | re.DOTALL,
    )
    analytical_processes = len(analytical_blocks)
    configured_thread_processes = sum(
        "cpus { moduleThreads(" in block
        or re.search(r"^\s+cpus 1$", block, re.MULTILINE) is not None
        for block in analytical_blocks
    )
    if configured_thread_processes != analytical_processes:
        errors.append(
            "Every analytical process must use its configured thread request"
        )
    configured_memory_processes = workflow_text.count("memory { moduleMemory(")
    if configured_memory_processes != analytical_processes:
        errors.append(
            "Every analytical process must use its configured memory request and command-line override"
        )

    entry_workflow = (ROOT / "nextflow" / "main_dsl2.nf").read_text()
    for entry in (
        "preprocessing",
        "assembly",
        "alignment",
        "coverage",
        "contact",
        "binning",
        "reassembly",
        "scaffolding",
        "annotation",
        "mge",
    ):
        if f"selected == '{entry}'" not in entry_workflow:
            errors.append(f"Missing descriptive Nextflow entry name: {entry}")

    launcher = (ROOT / "metahict").read_text()
    if not launcher.startswith("#!/bin/sh") or "conda info --base" not in launcher:
        errors.append("The primary launcher must bootstrap its Python runtime from Conda")
    if "metahict_manager.py" not in launcher or "exec" not in launcher:
        errors.append("The primary launcher must immediately delegate to the Python management interface")

    if errors:
        for error in errors:
            print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    print(
        "Architecture policy passed: Nextflow/Python own orchestration and stages; "
        "numbered module source trees are flat; shell is limited to bootstrap and packaging."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
