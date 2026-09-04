#!/usr/bin/env python3
"""User-facing management CLI for the METAHICT Nextflow workflow.

This module intentionally uses only the Python standard library so it works
before the project environments have been installed. Scientific stages are
Python entry points invoked directly by Nextflow.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from typing import Iterable, Optional, Sequence
from urllib.parse import unquote, urlparse


VERSION = "1.2.0"
PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "nextflow" / "assets" / "metahict_configuration.yaml"
LOCK_DIR = PROJECT_ROOT / "installation" / "locks" / "linux-64"
ENV_ROOT = PROJECT_ROOT / "conda_envs"
ENVIRONMENTS = (
    "metahict_nextflow_env",
    "metahict_env",
    "metabat2",
    "checkm2",
    "gtdbtk-2.4.0",
    "genomad",
    "ccfind_env",
)

ENTRY_MODULES = (
    "all",
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
)

REUSABLE_STAGE_INPUTS = (
    ("preprocessing_dir", "preprocessing_dir"),
    ("assembly_dir", "assembly_dir"),
    ("alignment_dir", "alignment_dir"),
    ("binning_dir", "binning_dir"),
    ("scaffolding_bin", "scaffolding_bin"),
    ("scaffolding_bam", "scaffolding_bam"),
    ("mge_fasta", "mge_fasta"),
    ("mge_alignment_dir", "mge_alignment_dir"),
    ("mge_contact_dir", "mge_contact_dir"),
)

DATABASES_BY_ENTRY = {
    "all": ("checkm", "checkm2", "gtdbtk", "genomad"),
    "binning": ("checkm", "checkm2"),
    "reassembly": ("checkm2",),
    "scaffolding": ("checkm2",),
    "annotation": ("gtdbtk",),
    "mge": ("genomad",),
}

DATABASE_ARGUMENTS = {
    "checkm": ("checkm_db", "checkm_db", "CheckM"),
    "checkm2": ("checkm2_db", "checkm2_db", "CheckM2"),
    "gtdbtk": ("gtdbtk_db", "gtdbtk_db", "GTDB-Tk"),
    "genomad": ("genomad_db", "genomad_db", "geNomad"),
}

ENZYME_ENTRY_MODULES = {"all", "contact", "binning", "scaffolding", "mge"}

LONG_READ_TYPES = (
    "pacbio-raw",
    "pacbio-corr",
    "pacbio-hifi",
    "nano-raw",
    "nano-corr",
    "nano-hq",
)

MODULE_HELP_DOCS = {
    "preprocessing": PROJECT_ROOT / "docs" / "modules" / "preprocessing.md",
    "assembly": PROJECT_ROOT / "docs" / "modules" / "assembly.md",
    "alignment": PROJECT_ROOT / "docs" / "modules" / "alignment.md",
    "coverage": PROJECT_ROOT / "docs" / "modules" / "coverage.md",
    "contact": PROJECT_ROOT / "docs" / "modules" / "contact.md",
    "binning": PROJECT_ROOT / "docs" / "modules" / "binning.md",
    "reassembly": PROJECT_ROOT / "docs" / "modules" / "reassembly.md",
    "scaffolding": PROJECT_ROOT / "docs" / "modules" / "scaffolding.md",
    "annotation": PROJECT_ROOT / "docs" / "modules" / "annotation.md",
    "mge": PROJECT_ROOT / "docs" / "modules" / "mge.md",
}

MODULE_RESOURCE_BEHAVIOR = {
    "preprocessing": (
        "Threads are passed to FastQC, BBDuk, and Clumpify. Nextflow passes "
        "80% of task memory to every BBTools Java command as -Xmx."
    ),
    "assembly": (
        "Threads are passed to the selected assembler and QUAST. Nextflow "
        "passes 80% of task memory to MEGAHIT or metaSPAdes; metaFlye has no "
        "native total-memory option, so its memory value is an executor allocation."
    ),
    "alignment": (
        "Threads are passed to BWA-MEM and SAMtools. BWA has no total-memory "
        "option, so task memory is an executor allocation; "
        "modules.alignment.sorting.memory_per_thread controls SAMtools memory per sort thread."
    ),
    "coverage": (
        "Threads are passed to BBMap and SAMtools. Nextflow passes 80% of task "
        "memory to BBMap as -Xmx."
    ),
    "contact": (
        "The current implementation is serial and therefore reserves one CPU. "
        "Its matrix operations have no native memory-limit option, so memory is "
        "an executor allocation."
    ),
    "binning": (
        "Threads are propagated through the binners, CheckM2, and refinement "
        "steps. These tools do not share one reliable total-memory flag, so task "
        "memory is an executor allocation."
    ),
    "reassembly": (
        "Threads are used by alignment, recruitment, CheckM2, and parallel "
        "SPAdes jobs. Nextflow passes 80% of task memory as the total SPAdes "
        "budget and divides it among concurrent bin assemblies."
    ),
    "scaffolding": (
        "Threads are passed to internal BWA/SAMtools alignment and CheckM2. "
        "YaHS has no total-memory flag, so task memory is an executor allocation; "
        "modules.scaffolding.alignment.sort_memory_per_thread controls SAMtools memory per thread."
    ),
    "annotation": (
        "Threads are passed to GTDB-Tk as --cpus and, unless overridden, as "
        "--pplacer_cpus. GTDB-Tk has no total-memory flag, so memory is an "
        "executor allocation."
    ),
    "mge": (
        "The shared thread setting is passed to geNomad, ccfind, and any generated "
        "Hi-C alignment; contact normalization remains serial. geNomad and ccfind "
        "have no shared total-memory flag, so memory is an executor allocation."
    ),
}


def required_databases_for_run(args: argparse.Namespace) -> tuple[str, ...]:
    entry_module = parse_entry_module(getattr(args, "entry_module", "all"))
    required = list(DATABASES_BY_ENTRY.get(entry_module, ()))
    if entry_module == "scaffolding" and getattr(
        args, "scaffolding_skip_checkm2", False
    ):
        required.remove("checkm2")
    return tuple(required)


class MetahictError(RuntimeError):
    """Expected user-facing failure."""


def is_macos_metadata_path(path: Path) -> bool:
    """Return whether a transferred path is macOS metadata rather than source."""
    return any(
        part == "__MACOSX" or part == ".DS_Store" or part.startswith("._")
        for part in path.parts
    )


class TeeTextIO(io.TextIOBase):
    """Write text to the user's terminal and an immutable run log."""

    def __init__(self, *streams: io.TextIOBase) -> None:
        self.streams = streams

    def write(self, text: str) -> int:
        for stream in self.streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()

    def isatty(self) -> bool:
        return any(stream.isatty() for stream in self.streams)


def display_command(command: Sequence[object]) -> str:
    return shlex.join(str(item) for item in command)


def run_command(
    command: Sequence[object],
    *,
    cwd: Optional[Path] = None,
    env: Optional[dict[str, str]] = None,
    capture: bool = False,
    announce: bool = True,
) -> subprocess.CompletedProcess[str]:
    argv = [str(item) for item in command]
    if announce:
        print(f"[RUN] {display_command(argv)}", flush=True)
    if not capture and (isinstance(sys.stdout, TeeTextIO) or isinstance(sys.stderr, TeeTextIO)):
        result = subprocess.run(
            argv,
            cwd=str(cwd) if cwd else None,
            env=env,
            check=False,
            text=True,
            capture_output=True,
        )
        if result.stdout:
            sys.stdout.write(result.stdout)
        if result.stderr:
            sys.stderr.write(result.stderr)
        if result.returncode:
            raise subprocess.CalledProcessError(
                result.returncode,
                argv,
                output=result.stdout,
                stderr=result.stderr,
            )
        return result
    return subprocess.run(
        argv,
        cwd=str(cwd) if cwd else None,
        env=env,
        check=True,
        text=True,
        capture_output=capture,
    )


def require_program(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise MetahictError(f"Required program was not found in PATH: {name}")
    return path


def explicit_urls(lines: Iterable[str]) -> list[str]:
    return sorted(
        line.strip()
        for line in lines
        if line.strip() and not line.lstrip().startswith("#") and line.strip() != "@EXPLICIT"
    )


def lock_urls(lock_file: Path) -> list[str]:
    if not lock_file.is_file():
        raise MetahictError(f"Missing environment lock: {lock_file}")
    return explicit_urls(lock_file.read_text().splitlines())


def parse_entry_module(value: str) -> str:
    """Return the canonical descriptive workflow entry name."""
    selected = value.strip().lower()
    if selected in ENTRY_MODULES:
        return selected
    valid = ", ".join(ENTRY_MODULES)
    raise argparse.ArgumentTypeError(f"unknown entry module '{value}'; choose one of: {valid}")


def positive_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def parse_memory_size(value: str) -> str:
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([MGT])(?:B)?\s*", value, re.IGNORECASE)
    if not match:
        raise argparse.ArgumentTypeError(
            "use a positive memory size such as 64000M, 64G, or 1T"
        )
    amount = float(match.group(1))
    if amount <= 0:
        raise argparse.ArgumentTypeError("memory must be greater than zero")
    number = str(int(amount)) if amount.is_integer() else str(amount)
    return f"{number} {match.group(2).upper()}B"


def _positive_integer_file(path: Path) -> int | None:
    try:
        value = path.read_text().strip()
    except (OSError, UnicodeError):
        return None
    if not value or value.lower() == "max":
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def local_cpu_capacity() -> int:
    """Return CPUs available to this local process, including common cgroup limits."""
    candidates: list[int] = []
    if hasattr(os, "sched_getaffinity"):
        try:
            candidates.append(len(os.sched_getaffinity(0)))
        except OSError:
            pass
    detected = os.cpu_count()
    if detected:
        candidates.append(detected)

    cpu_max = Path("/sys/fs/cgroup/cpu.max")
    try:
        quota_text, period_text = cpu_max.read_text().split()[:2]
        if quota_text != "max":
            quota, period = int(quota_text), int(period_text)
            if quota > 0 and period > 0:
                candidates.append(max(1, quota // period))
    except (OSError, UnicodeError, ValueError):
        pass
    quota = _positive_integer_file(
        Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
    )
    period = _positive_integer_file(
        Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    )
    if quota and period:
        candidates.append(max(1, quota // period))
    return max(1, min(candidates)) if candidates else 1


def local_memory_capacity_bytes() -> int:
    """Return physical RAM visible to the local process, capped by cgroups."""
    candidates: list[int] = []
    try:
        pages = int(os.sysconf("SC_PHYS_PAGES"))
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
        if pages > 0 and page_size > 0:
            candidates.append(pages * page_size)
    except (OSError, TypeError, ValueError):
        pass
    for path in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        limit = _positive_integer_file(path)
        if limit:
            candidates.append(limit)
    if not candidates:
        raise MetahictError("Could not determine memory available to the local executor")
    return min(candidates)


def local_resource_limits() -> tuple[int, int]:
    """Return local CPU and RAM ceilings used by Nextflow resourceLimits."""
    return local_cpu_capacity(), local_memory_capacity_bytes()


def _memory_bytes(value: str) -> int:
    match = re.fullmatch(
        r"\s*(\d+(?:\.\d+)?)\s*(B|K(?:I?B)?|M(?:I?B)?|G(?:I?B)?|T(?:I?B)?)\s*",
        str(value),
        re.IGNORECASE,
    )
    if not match:
        raise MetahictError(f"Invalid configured memory value: {value}")
    unit = match.group(2).upper()
    exponent = {
        "B": 0,
        "K": 1,
        "KB": 1,
        "KIB": 1,
        "M": 2,
        "MB": 2,
        "MIB": 2,
        "G": 3,
        "GB": 3,
        "GIB": 3,
        "T": 4,
        "TB": 4,
        "TIB": 4,
    }[unit]
    return int(float(match.group(1)) * (1024 ** exponent))


def _configured_resource_rows(config: Path) -> dict[str, dict[str, str]]:
    """Read the small resources map without adding a YAML dependency."""
    rows: dict[str, dict[str, str]] = {}
    in_resources = False
    current_module: str | None = None
    for raw_line in config.read_text().splitlines():
        content = raw_line.split(" #", 1)[0].rstrip()
        stripped = content.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(content) - len(stripped)
        if indent == 0:
            in_resources = stripped == "resources:"
            current_module = None
            continue
        if not in_resources:
            continue
        if indent == 2:
            module, value = stripped.split(":", 1)
            current_module = module.strip()
            rows.setdefault(current_module, {})
            inline = value.strip()
            if inline.startswith("{") and inline.endswith("}"):
                for item in inline[1:-1].split(","):
                    key, scalar = item.split(":", 1)
                    rows[current_module][key.strip()] = scalar.strip().strip("\"'")
            continue
        if indent == 4 and current_module and ":" in stripped:
            key, scalar = stripped.split(":", 1)
            rows[current_module][key.strip()] = scalar.strip().strip("\"'")
    return rows


def requested_resource_maxima(args: argparse.Namespace) -> tuple[int, int]:
    """Return the largest requested thread and memory values in this run."""
    defaults = _configured_resource_rows(DEFAULT_CONFIG)
    config = Path(getattr(args, "config", None) or DEFAULT_CONFIG).expanduser().resolve()
    configured = _configured_resource_rows(config)
    for module, values in configured.items():
        defaults.setdefault(module, {}).update(values)
    entry = parse_entry_module(getattr(args, "entry_module", "all"))
    selected = list(defaults) if entry == "all" else [entry]
    thread_override = getattr(args, "threads", None)
    memory_override = getattr(args, "memory", None)
    threads = [
        int(thread_override if thread_override is not None else defaults[module]["threads"])
        for module in selected
        if module in defaults
    ]
    memories = [
        _memory_bytes(
            memory_override if memory_override is not None else defaults[module]["memory"]
        )
        for module in selected
        if module in defaults
    ]
    return max(threads, default=1), max(memories, default=2 * 1024 ** 3)


def _format_capacity(memory_bytes: int) -> str:
    gibibytes = memory_bytes / (1024 ** 3)
    value = str(int(gibibytes)) if gibibytes.is_integer() else f"{gibibytes:.1f}"
    return f"{value} GB"


def warn_if_resources_are_capped(
    args: argparse.Namespace, cpu_limit: int, memory_limit: int
) -> None:
    requested_cpus, requested_memory = requested_resource_maxima(args)
    if requested_cpus <= cpu_limit and requested_memory <= memory_limit:
        return
    print(
        "[WARN] Requested resources exceed local capacity: "
        f"up to {requested_cpus} threads/{_format_capacity(requested_memory)} requested; "
        f"{cpu_limit} threads/{_format_capacity(memory_limit)} detected. "
        "Nextflow will cap each task to the detected limits and pass the effective "
        "values to its tools. Resource-intensive stages may run more slowly or need "
        "a larger machine.",
        file=sys.stderr,
    )


def configured_yaml_scalar(config: Path, *keys: str) -> str:
    """Read one scalar from the generated YAML without a Python dependency."""
    stack: list[tuple[int, str]] = []
    for raw_line in config.read_text().splitlines():
        stripped = raw_line.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(raw_line) - len(stripped)
        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        while stack and indent <= stack[-1][0]:
            stack.pop()
        path = tuple(item[1] for item in stack) + (key.strip(),)
        value = value.split(" #", 1)[0].strip()
        if not value:
            stack.append((indent, key.strip()))
            continue
        if path == keys:
            if value.lower() in {"", "null", "~"}:
                return ""
            if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
                try:
                    parsed = ast.literal_eval(value)
                except (SyntaxError, ValueError):
                    parsed = value[1:-1]
                return str(parsed).strip()
            return value
    return ""


def configuration_key_paths(config: Path) -> set[tuple[str, ...]]:
    """Return indentation-defined YAML mapping paths used by the configuration."""
    paths: set[tuple[str, ...]] = set()
    stack: list[tuple[int, str]] = []
    for raw_line in config.read_text().splitlines():
        stripped = raw_line.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        match = re.match(r"([A-Za-z_][A-Za-z0-9_-]*)\s*:", stripped)
        if not match:
            continue
        indent = len(raw_line) - len(stripped)
        while stack and indent <= stack[-1][0]:
            stack.pop()
        key = match.group(1)
        path = tuple(item[1] for item in stack) + (key,)
        paths.add(path)
        remainder = stripped[match.end():].split(" #", 1)[0].strip()
        if not remainder:
            stack.append((indent, key))
    return paths


def validate_configuration_schema(config: Path) -> None:
    """Reject obsolete or unknown keys that Nextflow would otherwise ignore."""
    section = ""
    obsolete: list[str] = []
    obsolete_module_keys = {
        "preprocessing": {
            "run_shotgun", "run_hic", "deduplicate_shotgun", "deduplicate_hic", "prefix",
            "minlen", "trimq", "qtrim", "ftl", "ftm", "ktrim", "k", "mink", "hdist",
            "adapter_ref", "skip_pre_qc_report", "skip_post_qc_report",
        },
        "assembly": {"min_len", "k_min", "k_max", "k_step", "merge_level", "k_list", "flye_method", "tmp_dir", "skip_quast", "keep_temp"},
        "alignment": {"samtools_filter", "bwa_options", "mapq", "min_intra_dist", "min_match_len", "sort_memory", "tmp_dir", "keep_sam", "skip_metrics"},
        "coverage": {"percent_identity", "min_mapq", "weight_mapq", "include_edge_bases", "max_edge_bases", "min_contig_length", "min_contig_depth", "bbmap_extra_args", "tmp_dir", "keep_temp"},
        "contact": {"metacc_min_signal", "metacc_min_len", "metacc_min_mapq", "metacc_min_match", "spurious_contact_percent", "coverage_file", "epsilon", "max_iter", "tol"},
        "binning": {"metacc_min_len", "metacc_min_signal", "metacc_min_mapq", "metacc_min_match", "metacc_min_binsize", "normcc_thres", "bin3c_min_len", "bin3c_min_signal", "bin3c_min_mapq", "bin3c_min_match", "bin3c_min_extent", "min_completeness", "max_contamination", "contamination_penalty", "min_input_bin_size", "max_input_bin_size", "binning_refiner_min_size", "tmp_dir", "keep_temp", "no_fasta", "no_report", "no_spades", "only_large", "skip_checkm2", "skip_refinement", "skip_consolidation", "keep_ambiguous", "remove_ambiguous", "seed"},
        "reassembly": {"cutoff_quantile", "top_k", "min_mapq", "min_match_len", "exclude_duplicates", "write_nonselected_hic", "min_contig_len", "strict_cut_off", "permissive_cut_off", "contamination_penalty", "skip_checkm2", "tmp_dir", "spades_mode", "spades_phred_offset", "spades_extra_args", "skip_residual_assembly", "keep_temp"},
        "scaffolding": {"resolution", "min_contig_len", "bwa_options", "samtools_filter", "sort_memory", "metacc_min_mapq", "metacc_min_len", "metacc_min_match", "metacc_min_signal", "yahs_resolutions", "yahs_min_mapq", "yahs_min_contig_len", "yahs_rounds", "yahs_no_contig_ec", "yahs_no_scaffold_ec", "yahs_no_mem_check", "yahs_extra_args", "normcc_thres", "heatmap_max_image", "skip_checkm2", "tmp_dir", "keep_temp"},
        "annotation": {"pplacer_cpus", "extension", "prefix", "skip_ani_screen", "mash_db", "min_perc_aa", "min_af", "full_tree", "scratch_dir", "tmp_dir", "force", "gtdbtk_extra_args"},
        "mge": {"genomad_splits", "genomad_sensitivity", "genomad_cleanup", "genomad_restart", "genomad_preset", "genomad_min_score", "genomad_max_fdr", "genomad_extra_args", "association_filter", "zscore_threshold", "fixed_contact_threshold", "top_percent", "min_raw_contacts", "ccfind_terminal_fragment_size", "ccfind_min_identity", "ccfind_min_aligned_length", "min_contact_strength", "tmp_dir"},
    }
    nested_obsolete_paths = {
        ("modules", "mge", "alignment", key)
        for key in obsolete_module_keys["alignment"]
    } | {
        ("modules", "mge", "contact", key)
        for key in obsolete_module_keys["contact"]
    }
    stack: list[tuple[int, str]] = []
    for line_number, raw_line in enumerate(config.read_text().splitlines(), start=1):
        stripped = raw_line.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(raw_line) - len(stripped)
        if ":" in stripped:
            key, value = stripped.split(":", 1)
            while stack and indent <= stack[-1][0]:
                stack.pop()
            path = tuple(item[1] for item in stack) + (key.strip(),)
            if len(path) == 3 and path[0] == "modules":
                module_name, module_key = path[1], path[2]
                if module_key in obsolete_module_keys.get(module_name, set()):
                    obsolete.append(f"line {line_number}: {'.'.join(path)}")
            if path in nested_obsolete_paths:
                obsolete.append(f"line {line_number}: {'.'.join(path)}")
            if not value.strip():
                stack.append((indent, key.strip()))
        if indent == 0 and ":" in stripped:
            section = stripped.split(":", 1)[0].strip()
            if section == "workflow":
                obsolete.append(f"line {line_number}: workflow")
            continue
        if section == "resources" and re.search(
            r"(?:^|[,{]\s*)(cpus|time)\s*:", stripped
        ):
            key = re.search(r"(?:^|[,{]\s*)(cpus|time)\s*:", stripped).group(1)
            obsolete.append(f"line {line_number}: resources ... {key}")
        if (
            section == "resources"
            and indent == 2
            and stripped.split(":", 1)[0] in {"mge_alignment", "mge_contact"}
        ):
            key = stripped.split(":", 1)[0]
            obsolete.append(f"line {line_number}: resources.{key}")
        if section == "modules" and re.match(
            r"preprocessing_(?:sg|hic)\s*:", stripped
        ):
            key = stripped.split(":", 1)[0]
            obsolete.append(f"line {line_number}: modules.{key}")
        if (
            section == "modules"
            and indent == 2
            and stripped.split(":", 1)[0] in {"mge_alignment", "mge_contact"}
        ):
            key = stripped.split(":", 1)[0]
            obsolete.append(f"line {line_number}: modules.{key}")
    template_paths = configuration_key_paths(DEFAULT_CONFIG)
    configured_paths = configuration_key_paths(config)
    resource_modules = {
        path[1]
        for path in template_paths
        if len(path) == 2 and path[0] == "resources"
    }
    allowed_resource_leaves = {
        ("resources", module, key)
        for module in resource_modules
        for key in ("threads", "memory")
    }
    unknown = sorted(
        path
        for path in configured_paths
        if path not in template_paths and path not in allowed_resource_leaves
    )
    obsolete.extend(f"unknown key: {'.'.join(path)}" for path in unknown)
    if obsolete:
        raise MetahictError(
            "Obsolete configuration schema detected:\n  "
            + "\n  ".join(obsolete)
            + "\nCreate the current schema with './metahict config --force' after "
            "backing up your existing metahict_configuration.yaml."
        )


def environment_prefix(name: str, *, allow_legacy: bool = True) -> Path:
    """Resolve one environment, including the legacy Java-environment name."""
    prefix = ENV_ROOT / name
    if prefix.is_dir() or not allow_legacy or name != "metahict_nextflow_env":
        return prefix
    legacy = ENV_ROOT / "metahict_venv"
    if legacy.is_dir():
        return legacy
    return prefix


def verify_lock_checksums(*, verbose: bool = True) -> int:
    checksum_file = LOCK_DIR.parent / "SHA256SUMS"
    if not checksum_file.is_file():
        raise MetahictError(f"Missing checksum manifest: {checksum_file}")
    verified = 0
    for line in checksum_file.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = checksum_file.parent / relative.lstrip("* ")
        observed = hashlib.sha256(target.read_bytes()).hexdigest()
        if observed != expected:
            raise MetahictError(f"Lock checksum mismatch: {target}")
        verified += 1
        if verbose:
            print(f"[PASS] lock checksum: {target.name}")
    return verified


def installed_urls(env_prefix: Path, *, verbose: bool = True) -> list[str]:
    result = run_command(
        ["conda", "list", "--explicit", "-p", env_prefix],
        capture=True,
        announce=verbose,
    )
    return explicit_urls(result.stdout.splitlines())


def verify_environment(
    name: str, *, allow_legacy: bool = True, verbose: bool = True
) -> None:
    prefix = environment_prefix(name, allow_legacy=allow_legacy)
    if not prefix.is_dir():
        raise MetahictError(f"Missing environment: {prefix}")
    expected = lock_urls(LOCK_DIR / f"{name}.explicit.txt")
    observed = installed_urls(prefix, verbose=verbose)
    if observed != expected:
        raise MetahictError(
            f"Environment '{name}' does not match its exact package lock. "
            f"Recreate {prefix}."
        )
    if prefix.name != name:
        print(
            f"[COMPATIBILITY] using legacy environment {prefix}; "
            "run './metahict install' to create conda_envs/metahict_nextflow_env"
        )
    if verbose:
        print(f"[PASS] environment lock: {name} ({prefix})")


def verify_runtime(*, verbose: bool = True) -> None:
    if platform.system() != "Linux":
        raise MetahictError("The distributed environment locks support Linux only.")
    require_program("conda")
    checksum_count = verify_lock_checksums(verbose=verbose)
    for name in ENVIRONMENTS:
        verify_environment(name, verbose=verbose)

    run_command(
        [
            "conda",
            "run",
            "-p",
            ENV_ROOT / "metahict_env",
            "python",
            "-c",
            (
                "from importlib.metadata import version; "
                "assert version('Binning-refiner') == '1.4.3'"
            ),
        ],
        announce=verbose,
    )
    summarizer = ENV_ROOT / "metabat2" / "bin" / "jgi_summarize_bam_contig_depths"
    if not os.access(summarizer, os.X_OK):
        raise MetahictError(f"MetaBAT2 coverage executable is missing: {summarizer}")
    for name in ("checkm2", "diamond"):
        executable = ENV_ROOT / "checkm2" / "bin" / name
        if not os.access(executable, os.X_OK):
            raise MetahictError(f"CheckM2 executable is missing: {executable}")
    adapter_reference = ENV_ROOT / "metahict_env" / "share" / "bbmap" / "resources" / "adapters.fa"
    if not adapter_reference.is_file():
        raise MetahictError(f"BBTools adapter reference is missing: {adapter_reference}")
    ccfind = PROJECT_ROOT / "external" / "bin" / "ccfind"
    if not os.access(ccfind, os.X_OK):
        raise MetahictError(f"ccfind launcher is missing: {ccfind}")
    ssearch36 = ENV_ROOT / "ccfind_env" / "bin" / "ssearch36"
    if not os.access(ssearch36, os.X_OK):
        raise MetahictError(f"ccfind companion executable is missing: {ssearch36}")
    if verbose:
        print(
            "[PASS] pinned Binning_refiner, MetaBAT2, CheckM2/DIAMOND, "
            "BBTools resources, and ccfind/FASTA-suite tools"
        )
    else:
        print(
            "[PASS] Runtime preflight: "
            f"{checksum_count} lock checksums, {len(ENVIRONMENTS)} environments, "
            "and pinned resources verified"
        )


def create_environment(name: str) -> None:
    lock_file = LOCK_DIR / f"{name}.explicit.txt"
    prefix = ENV_ROOT / name
    if prefix.exists():
        verify_environment(name, allow_legacy=False)
        return
    prefix.parent.mkdir(parents=True, exist_ok=True)
    run_command(["conda", "create", "-y", "-p", prefix, "--file", lock_file])
    verify_environment(name, allow_legacy=False)


def download(url: str, destination: Path) -> None:
    require_program("curl")
    destination.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "curl",
            "--fail",
            "--location",
            "--continue-at",
            "-",
            "--output",
            destination,
            url,
        ]
    )


def extract_tar_gz(archive: Path, destination: Path) -> None:
    require_program("tar")
    destination.mkdir(parents=True, exist_ok=True)
    run_command(["tar", "-xzf", archive, "-C", destination])


def install_ccfind() -> None:
    version = "1.4.7"
    commit = "674366b49dd31cb909c2e52834e4ec8ede8919e7"
    expected_sha = "abee22a03ab2e7c475d08361f90fa409a5626733649e4eb39b0cdbc7036ff463"
    external = PROJECT_ROOT / "external"
    archive = external / f"ccfind-{version}.tar.gz"
    source = external / f"ccfind-{commit}"
    binary_dir = external / "bin"
    url = f"https://codeload.github.com/yosuken/ccfind/tar.gz/{commit}"

    if not archive.is_file():
        download(url, archive)
    observed_sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    if observed_sha != expected_sha:
        raise MetahictError(f"ccfind checksum mismatch: {archive}")
    if not source.is_dir():
        extract_tar_gz(archive, external)

    binary = source / "ccfind"
    binary.chmod(binary.stat().st_mode | 0o111)
    binary_dir.mkdir(parents=True, exist_ok=True)
    launcher = binary_dir / "ccfind"
    ccfind_environment_bin = ENV_ROOT / "ccfind_env" / "bin"
    launcher.write_text(
        "#!/usr/bin/env python3\n"
        "import os, sys\n"
        f"os.environ['PATH'] = {str(ccfind_environment_bin)!r} + os.pathsep + os.environ.get('PATH', '')\n"
        f"os.execv({str(binary)!r}, [{str(binary)!r}, *sys.argv[1:]])\n"
    )
    launcher.chmod(0o755)
    link = external / "ccfind"
    if link.is_symlink() or link.exists():
        if link.is_symlink() and link.resolve() == source.resolve():
            pass
        else:
            raise MetahictError(f"Refusing to replace existing path: {link}")
    else:
        link.symlink_to(source, target_is_directory=True)
    print(f"[PASS] ccfind {version}: {binary}")


def command_install(_: argparse.Namespace) -> None:
    if platform.system() != "Linux":
        raise MetahictError("Installation is supported only on Linux for this release.")
    require_program("conda")
    require_program("git")
    verify_lock_checksums()
    for name in ENVIRONMENTS:
        create_environment(name)
    install_ccfind()
    run_command(
        [
            "conda",
            "run",
            "-p",
            ENV_ROOT / "metahict_env",
            "python",
            "-m",
            "pip",
            "install",
            "--no-deps",
            "--no-build-isolation",
            "--upgrade",
            "--force-reinstall",
            "-r",
            PROJECT_ROOT / "installation" / "pip-requirements.txt",
        ]
    )
    verify_runtime()
    print("[PASS] METAHICT locked runtime installation")


def database_directory(target: str, path: Optional[str]) -> Path:
    if path:
        return Path(path).expanduser().resolve()
    defaults = {
        "checkm": PROJECT_ROOT / "databases" / "checkm_db",
        "checkm2": PROJECT_ROOT / "databases" / "checkm2_db",
        "gtdbtk": PROJECT_ROOT / "databases" / "gtdbtk_db",
        "genomad": PROJECT_ROOT / "databases" / "genomad_db",
    }
    return defaults[target]


def install_checkm_database(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination / "checkm_data.tar.gz"
    download(
        "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz",
        archive,
    )
    extract_tar_gz(archive, destination)
    run_command(
        [
            "conda",
            "run",
            "-p",
            ENV_ROOT / "metahict_env",
            "checkm",
            "data",
            "setRoot",
            destination,
        ]
    )
    archive.unlink(missing_ok=True)


def install_checkm2_database(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "conda",
            "run",
            "-p",
            ENV_ROOT / "checkm2",
            "checkm2",
            "database",
            "--download",
            "--path",
            destination,
        ]
    )
    expected = destination / "CheckM2_database" / "uniref100.KO.1.dmnd"
    if not expected.is_file() or expected.stat().st_size == 0:
        raise MetahictError(f"Incomplete CheckM2 database: {expected}")


def install_gtdbtk_database(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination / "gtdbtk_data.tar.gz"
    download(
        "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/"
        "auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz",
        archive,
    )
    extract_tar_gz(archive, destination)
    release = destination / "release220"
    if not release.is_dir():
        raise MetahictError(f"Incomplete GTDB-Tk database: {release}")
    run_command(
        [
            "conda",
            "env",
            "config",
            "vars",
            "set",
            "-p",
            ENV_ROOT / "gtdbtk-2.4.0",
            f"GTDBTK_DATA_PATH={release}",
        ]
    )
    archive.unlink(missing_ok=True)


def install_genomad_database(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "conda",
            "run",
            "-p",
            ENV_ROOT / "genomad",
            "genomad",
            "download-database",
            destination,
        ]
    )
    expected = destination / "genomad_db" / "version.txt"
    if not expected.is_file():
        raise MetahictError(f"Incomplete geNomad database: {expected}")


def command_database(args: argparse.Namespace) -> None:
    require_program("conda")
    targets = ("checkm", "checkm2", "gtdbtk", "genomad") if args.target == "all" else (args.target,)
    if args.target == "all" and args.path:
        raise MetahictError("--path is valid for one database; use --root with 'all'.")
    installers = {
        "checkm": install_checkm_database,
        "checkm2": install_checkm2_database,
        "gtdbtk": install_gtdbtk_database,
        "genomad": install_genomad_database,
    }
    for target in targets:
        destination = (
            Path(args.root).expanduser().resolve() / database_directory(target, None).name
            if args.target == "all" and args.root
            else database_directory(target, args.path)
        )
        env_name = {
            "checkm": "metahict_env",
            "checkm2": "checkm2",
            "gtdbtk": "gtdbtk-2.4.0",
            "genomad": "genomad",
        }[target]
        if not (ENV_ROOT / env_name).is_dir():
            raise MetahictError("Install the locked runtime before downloading databases.")
        print(f"[INFO] Installing {target} database to {destination}")
        installers[target](destination)
        print(f"[PASS] {target} database: {destination}")


def database_status() -> list[tuple[str, bool, Path]]:
    checks = [
        ("CheckM", (PROJECT_ROOT / "databases" / "checkm_db").is_dir(), PROJECT_ROOT / "databases" / "checkm_db"),
        (
            "CheckM2",
            (PROJECT_ROOT / "databases" / "checkm2_db" / "CheckM2_database" / "uniref100.KO.1.dmnd").is_file(),
            PROJECT_ROOT / "databases" / "checkm2_db" / "CheckM2_database" / "uniref100.KO.1.dmnd",
        ),
        ("GTDB-Tk", (PROJECT_ROOT / "databases" / "gtdbtk_db" / "release220").is_dir(), PROJECT_ROOT / "databases" / "gtdbtk_db" / "release220"),
        ("geNomad", (PROJECT_ROOT / "databases" / "genomad_db" / "genomad_db" / "version.txt").is_file(), PROJECT_ROOT / "databases" / "genomad_db" / "genomad_db"),
    ]
    return checks


def command_doctor(args: argparse.Namespace) -> None:
    failures = 0

    def report(label: str, passed: bool, detail: object) -> None:
        nonlocal failures
        print(f"[{'PASS' if passed else 'FAIL'}] {label}: {detail}")
        if not passed:
            failures += 1

    report("operating system", platform.system() == "Linux", platform.platform())
    report("architecture", platform.machine() in {"x86_64", "AMD64"}, platform.machine())
    report(
        "Conda-provided Python runtime",
        sys.version_info >= (3, 9),
        platform.python_version(),
    )
    for program in ("conda", "curl", "tar", "git"):
        report(program, shutil.which(program) is not None, shutil.which(program) or "not found")

    try:
        verify_lock_checksums()
    except (MetahictError, OSError) as error:
        report("environment lock checksums", False, error)
    else:
        report("environment lock checksums", True, LOCK_DIR.parent / "SHA256SUMS")

    if args.runtime:
        for name in ENVIRONMENTS:
            try:
                verify_environment(name)
            except (MetahictError, subprocess.CalledProcessError) as error:
                report(f"environment {name}", False, error)
            else:
                report(f"environment {name}", True, environment_prefix(name))

    if args.databases:
        for label, passed, path in database_status():
            report(f"database {label}", passed, path)

    if failures:
        raise MetahictError(f"Doctor found {failures} failed check(s).")


def source_checks() -> None:
    run_command([sys.executable, PROJECT_ROOT / "nextflow" / "ci" / "check_no_vendored_tools.py"])
    run_command([sys.executable, PROJECT_ROOT / "nextflow" / "ci" / "check_architecture_policy.py"])
    run_command([sys.executable, PROJECT_ROOT / "nextflow" / "ci" / "check_documentation.py"])
    run_command([sys.executable, PROJECT_ROOT / "nextflow" / "ci" / "check_example_dataset.py"])
    run_command(
        [
            sys.executable,
            "-m",
            "unittest",
            str(PROJECT_ROOT / "nextflow" / "ci" / "test_binning_refiner_adapter.py"),
        ]
    )
    run_command(
        [
            sys.executable,
            "-m",
            "unittest",
            str(PROJECT_ROOT / "nextflow" / "ci" / "test_manager_cli.py"),
        ]
    )
    run_command(
        [
            sys.executable,
            "-m",
            "unittest",
            str(PROJECT_ROOT / "nextflow" / "ci" / "test_stage_interfaces.py"),
        ]
    )
    run_command(
        [
            sys.executable,
            "-m",
            "unittest",
            str(PROJECT_ROOT / "nextflow" / "ci" / "test_scientific_comparison.py"),
        ]
    )
    for path in PROJECT_ROOT.rglob("*.py"):
        if is_macos_metadata_path(path) or any(
            part in {".git", "conda_envs", "nextflow_work", "work"}
            for part in path.parts
        ):
            continue
        ast.parse(path.read_text(), filename=str(path))
    print("[PASS] Python syntax")
    run_command([require_program("sh"), "-n", PROJECT_ROOT / "metahict"])
    print("[PASS] Conda/Python bootstrap syntax")
    shell_files = sorted(
        path
        for path in PROJECT_ROOT.rglob("*.sh")
        if not is_macos_metadata_path(path)
        and not any(
            part in {".git", "conda_envs", "nextflow_work", "work"}
            for part in path.parts
        )
    )
    if shell_files:
        run_command([require_program("bash"), "-n", *shell_files])
    print("[PASS] shell bootstrap/packaging syntax")
    verify_lock_checksums()


def nextflow_environment() -> dict[str, str]:
    env = os.environ.copy()
    nextflow_prefix = environment_prefix("metahict_nextflow_env")
    nextflow_bin = nextflow_prefix / "bin"
    env["NXF_HOME"] = str(PROJECT_ROOT / "nextflow" / ".nextflow")
    path_entries = [str(PROJECT_ROOT / "nextflow" / "bin")]
    if (nextflow_bin / "java").is_file():
        env["JAVA_HOME"] = str(nextflow_prefix)
        path_entries.extend(
            [
                str(nextflow_bin),
                str(ENV_ROOT / "metahict_env" / "bin"),
                str(PROJECT_ROOT / "external" / "bin"),
            ]
        )
    path_entries.append(env.get("PATH", ""))
    env["PATH"] = os.pathsep.join(path_entries)
    return env


def run_stub_test(profile: str) -> None:
    if profile == "local":
        verify_runtime()
    nextflow = PROJECT_ROOT / "nextflow" / "bin" / "nextflow"
    if not nextflow.is_file():
        raise MetahictError(f"Bundled Nextflow launcher is missing: {nextflow}")
    with tempfile.TemporaryDirectory(prefix=f"metahict-{profile}-smoke-") as temporary:
        root = Path(temporary)
        dummy = root / "dummy_db"
        for directory in ("checkm_db", "gtdbtk_db", "genomad_db"):
            (dummy / directory).mkdir(parents=True, exist_ok=True)
        (dummy / "checkm2.dmnd").touch()
        stub_input_directory = root / "stub_inputs"
        run_command(
            [
                sys.executable,
                PROJECT_ROOT / "nextflow" / "ci" / "create_stub_inputs.py",
                "--output-dir",
                stub_input_directory,
            ]
        )
        stub_samplesheet = stub_input_directory / "samplesheet.csv"
        results = root / "results"
        command = [
            nextflow,
            "run",
            PROJECT_ROOT / "nextflow" / "main_dsl2.nf",
            "-params-file",
            PROJECT_ROOT / "nextflow" / "assets" / "example_dataset_configuration.yaml",
            "-profile",
            profile,
            "-c",
            PROJECT_ROOT / "nextflow" / "ci" / "stub_resources.config",
            "--samplesheet",
            stub_samplesheet,
            "--out_root",
            results,
            "--report_dir",
            root / "reports",
            "-work-dir",
            root / "work",
            "--checkm_db",
            dummy / "checkm_db",
            "--checkm2_db",
            dummy / "checkm2.dmnd",
            "--gtdbtk_db",
            dummy / "gtdbtk_db",
            "--genomad_db",
            dummy / "genomad_db",
            "--threads",
            "2",
            "--clean",
            "true",
            "--chain",
            "true",
            "-stub-run",
            "-ansi-log",
            "false",
        ]
        run_command(command, cwd=PROJECT_ROOT, env=nextflow_environment())
        unexpected_scaffolding = sorted(results.glob("*/8_scaffolding"))
        if unexpected_scaffolding:
            raise MetahictError(
                "The default complete workflow unexpectedly ran scaffolding:\n  "
                + "\n  ".join(str(path) for path in unexpected_scaffolding)
            )
        run_stub_scaffolding_entries(
            profile=profile,
            source_samplesheet=stub_samplesheet,
            results=results,
            root=root,
            nextflow=nextflow,
            checkm2_database=dummy / "checkm2.dmnd",
        )
        run_command(
            [
                sys.executable,
                PROJECT_ROOT / "nextflow" / "bin" / "check_expected_outputs.py",
                "--root",
                results,
                "--manifest",
                PROJECT_ROOT / "nextflow" / "tests" / "expected" / "workflow_stub_outputs.tsv",
            ]
        )
        unexpected_long_read_outputs = (
            results / "long_read_example" / "1_preprocessing" / "sg",
            results / "long_read_example" / "7_reassembly",
        )
        present = [str(path) for path in unexpected_long_read_outputs if path.exists()]
        if present:
            raise MetahictError(
                "Long-read stub unexpectedly ran a short-read-only stage:\n  "
                + "\n  ".join(present)
            )
    print(f"[PASS] Nextflow {profile} stub test")


def sample_rows(samplesheet: Path) -> list[dict[str, str]]:
    """Read non-empty sample rows from a METAHICT samplesheet."""
    with samplesheet.open(newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle)
            if (row.get("sample") or "").strip()
        ]


def write_one_sample_sheet(
    source_samplesheet: Path, destination: Path, selected_row: dict[str, str]
) -> None:
    """Write one existing samplesheet row without changing its columns."""
    with source_samplesheet.open(newline="") as source:
        fieldnames = csv.DictReader(source).fieldnames
    if not fieldnames:
        raise MetahictError(f"Samplesheet has no header: {source_samplesheet}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({field: selected_row.get(field, "") for field in fieldnames})


def example_scaffolding_bins(results: Path, row: dict[str, str]) -> list[Path]:
    """Return all final MAG FASTAs available for optional example scaffolding."""
    sample = row["sample"].strip()
    if (row.get("long_read_type") or "").strip():
        directory = results / sample / "6_binning" / "metahict" / "final_bins"
    else:
        directory = results / sample / "7_reassembly" / "reassembled_bins"
    if not directory.is_dir():
        return []
    suffixes = {".fa", ".fasta", ".fna"}
    return sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in suffixes
    )


def run_stub_scaffolding_entries(
    *,
    profile: str,
    source_samplesheet: Path,
    results: Path,
    root: Path,
    nextflow: Path,
    checkm2_database: Path,
) -> None:
    """Exercise the standalone scaffolding entry after the default stub run."""
    for row in sample_rows(source_samplesheet):
        sample = row["sample"].strip()
        bins = example_scaffolding_bins(results, row)
        if not bins:
            print(f"[WARN] Stub sample {sample} produced no bin for scaffolding")
            continue
        one_sample_sheet = root / "scaffolding_samplesheets" / f"{sample}.csv"
        write_one_sample_sheet(source_samplesheet, one_sample_sheet, row)
        for bin_fasta in bins:
            safe_bin = safe_log_name(bin_fasta.stem)
            command = [
                nextflow,
                "run",
                PROJECT_ROOT / "nextflow" / "main_dsl2.nf",
                "-params-file",
                PROJECT_ROOT / "nextflow" / "assets" / "example_dataset_configuration.yaml",
                "-profile",
                profile,
                "-c",
                PROJECT_ROOT / "nextflow" / "ci" / "stub_resources.config",
                "--entry_module",
                "scaffolding",
                "--samplesheet",
                one_sample_sheet,
                "--scaffolding_bin",
                bin_fasta,
                "--preprocessing_dir",
                results / sample / "1_preprocessing",
                "--out_root",
                results,
                "--report_dir",
                root / "scaffolding_reports" / sample / safe_bin,
                "-work-dir",
                root / "scaffolding_work" / sample / safe_bin,
                "--checkm2_db",
                checkm2_database,
                "--threads",
                "2",
                "-stub-run",
                "-ansi-log",
                "false",
            ]
            run_command(command, cwd=PROJECT_ROOT, env=nextflow_environment())


def validate_example_scaffolding(
    results: Path, attempts: Sequence[tuple[str, Path]]
) -> None:
    """Validate interfaces without requiring a particular biological outcome."""
    for sample, bin_fasta in attempts:
        output = results / sample / "8_scaffolding" / bin_fasta.stem
        status_file = output / "scaffolding_status.tsv"
        if not status_file.is_file() or status_file.stat().st_size == 0:
            raise MetahictError(
                f"Scaffolding did not report a status for {sample}:{bin_fasta.name}"
            )
        with status_file.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(rows) != 1 or rows[0].get("status") not in {"completed", "skipped"}:
            raise MetahictError(f"Invalid scaffolding status file: {status_file}")


def append_example_database_options(command: list[str], args: argparse.Namespace) -> None:
    for argument_name, _, _ in DATABASE_ARGUMENTS.values():
        value = getattr(args, argument_name, None)
        if value:
            command.extend(
                [f"--{argument_name.replace('_', '-')}", str(Path(value).expanduser())]
            )


def run_example_test(args: argparse.Namespace) -> None:
    """Run the bundled data through core stages and optional scaffolding."""
    results = Path(args.outdir).expanduser().resolve()
    samplesheet = PROJECT_ROOT / "nextflow" / "assets" / "example_dataset_samplesheet.csv"
    config = Path(args.config).expanduser().resolve()
    base = [
        "run",
        "--samplesheet",
        str(samplesheet),
        "--config",
        str(config),
        "--outdir",
        str(results),
    ]
    if args.resume:
        base.append("--resume")
    append_example_database_options(base, args)
    command_run(build_parser().parse_args(base))

    attempts: list[tuple[str, Path]] = []
    for row in sample_rows(samplesheet):
        sample = row["sample"].strip()
        bins = example_scaffolding_bins(results, row)
        if not bins:
            print(f"[WARN] No final MAGs were available for scaffolding in {sample}")
            continue
        one_sample_sheet = (
            results
            / "nextflow_reports"
            / "example_test_samplesheets"
            / f"{sample}.csv"
        )
        write_one_sample_sheet(samplesheet, one_sample_sheet, row)
        for bin_fasta in bins:
            scaffold = [
                "run",
                "--entry-module",
                "scaffolding",
                "--samplesheet",
                str(one_sample_sheet),
                "--config",
                str(config),
                "--scaffolding-bin",
                str(bin_fasta),
                "--preprocessing-dir",
                str(results / sample / "1_preprocessing"),
                "--outdir",
                str(results),
            ]
            if args.resume:
                scaffold.append("--resume")
            append_example_database_options(scaffold, args)
            command_run(build_parser().parse_args(scaffold))
            attempts.append((sample, bin_fasta))

    run_command(
        [
            sys.executable,
            PROJECT_ROOT / "nextflow" / "bin" / "check_expected_outputs.py",
            "--root",
            results,
            "--manifest",
            PROJECT_ROOT / "nextflow" / "tests" / "expected" / "example_dataset_outputs.tsv",
        ]
    )
    validate_example_scaffolding(results, attempts)
    print("[PASS] Bundled example test completed")


def command_test(args: argparse.Namespace) -> None:
    if args.scope in {"source", "all"}:
        source_checks()
    if args.scope in {"workflow", "all"}:
        run_stub_test("stub")
    if args.scope == "example":
        run_example_test(args)


def command_compare(args: argparse.Namespace) -> None:
    baseline = Path(args.baseline).expanduser().resolve()
    candidate = Path(args.candidate).expanduser().resolve()
    manifest = (
        Path(args.manifest).expanduser().resolve()
        if args.manifest
        else PROJECT_ROOT / "nextflow" / "tests" / "scientific_regression_manifest.tsv"
    )
    report = Path(args.report).expanduser().resolve()
    for label, path in (("baseline", baseline), ("candidate", candidate)):
        if not path.is_dir():
            raise MetahictError(f"Scientific comparison {label} directory not found: {path}")
    if not manifest.is_file():
        raise MetahictError(f"Scientific comparison manifest not found: {manifest}")
    run_command(
        [
            sys.executable,
            PROJECT_ROOT / "nextflow" / "bin" / "compare_scientific_outputs.py",
            "--baseline",
            baseline,
            "--candidate",
            candidate,
            "--manifest",
            manifest,
            "--report",
            report,
        ]
    )


def default_database_paths() -> dict[str, Path]:
    return {
        "checkm": PROJECT_ROOT / "databases" / "checkm_db",
        "checkm2": PROJECT_ROOT / "databases" / "checkm2_db" / "CheckM2_database" / "uniref100.KO.1.dmnd",
        "gtdbtk": PROJECT_ROOT / "databases" / "gtdbtk_db" / "release220",
        "genomad": PROJECT_ROOT / "databases" / "genomad_db" / "genomad_db",
    }


REPORT_LINKS = (
    "run.log",
    "parameters.json",
    "provenance.json",
    "software_versions.tsv",
    "environment_locks.tsv",
    "samplesheet.sha256",
    "configuration.sha256",
    "trace.txt",
    "report.html",
    "timeline.html",
    "dag.html",
    "nextflow.log",
    "process_logs",
    "failure_summary.txt",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def report_root(args: argparse.Namespace) -> Path:
    out_root = Path(args.outdir).expanduser().resolve()
    return Path(args.report_dir).expanduser().resolve() if args.report_dir else out_root / "nextflow_reports"


def work_root(args: argparse.Namespace) -> Path:
    out_root = Path(args.outdir).expanduser().resolve()
    return Path(args.work_dir).expanduser().resolve() if args.work_dir else out_root / "nextflow_work"


def create_run_directory(reports: Path, now: datetime | None = None) -> tuple[str, Path]:
    timestamp = (now or datetime.now(timezone.utc)).astimezone(timezone.utc)
    base = timestamp.strftime("%Y%m%dT%H%M%SZ")
    runs = reports / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    for suffix in range(1000):
        run_id = base if suffix == 0 else f"{base}-{suffix:03d}"
        run_dir = runs / run_id
        try:
            run_dir.mkdir()
        except FileExistsError:
            continue
        (run_dir / "process_logs").mkdir()
        return run_id, run_dir
    raise MetahictError(f"Could not allocate a unique run directory below {runs}")


def replace_report_link(link: Path, target: Path, legacy_dir: Path) -> None:
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        legacy_dir.mkdir(parents=True, exist_ok=True)
        legacy_target = legacy_dir / link.name
        counter = 1
        while legacy_target.exists():
            legacy_target = legacy_dir / f"{link.name}.{counter}"
            counter += 1
        shutil.move(str(link), str(legacy_target))
    relative_target = os.path.relpath(target, start=link.parent)
    link.symlink_to(relative_target, target_is_directory=target.is_dir())


def prepare_report_links(reports: Path, run_id: str, run_dir: Path) -> str | None:
    previous_run: str | None = None
    latest = reports / "latest"
    if latest.is_symlink():
        try:
            previous_run = latest.resolve(strict=True).name
        except FileNotFoundError:
            previous_run = None
    legacy_dir = reports / "runs" / f"legacy-before-{run_id}"
    replace_report_link(latest, run_dir, legacy_dir)
    for name in REPORT_LINKS:
        replace_report_link(reports / name, run_dir / name, legacy_dir)
    return previous_run


def atomic_write_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    os.replace(temporary, path)


def git_revision() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    return result.stdout.strip() or None


def bundled_nextflow_version() -> str:
    launcher = PROJECT_ROOT / "nextflow" / "bin" / "nextflow"
    match = re.search(r"^NXF_VER=.*?'([^']+)'", launcher.read_text(), flags=re.MULTILINE)
    return match.group(1) if match else "unknown"


def locked_package(url: str) -> tuple[str, str, str] | None:
    filename = unquote(Path(urlparse(url).path).name)
    for suffix in (".tar.bz2", ".conda"):
        if filename.endswith(suffix):
            filename = filename[: -len(suffix)]
            break
    fields = filename.rsplit("-", 2)
    return tuple(fields) if len(fields) == 3 else None


def pinned_python_requirement(line: str) -> tuple[str, str, str] | None:
    if " @ " not in line:
        return None
    name, source = (field.strip() for field in line.split(" @ ", 1))
    if source.startswith("git+") and "@" in source:
        return name, source.rsplit("@", 1)[1], "git-commit"
    wheel = re.search(r"[-_]([0-9][A-Za-z0-9.]*)-[^-/#]+-[^-/#]+-[^/#]+\.whl", source)
    if wheel:
        return name, wheel.group(1), "wheel"
    return name, "pinned", "direct-reference"


def write_environment_metadata(run_dir: Path) -> None:
    lock_lines = ["environment\tprefix\tlock_file\tsha256"]
    version_lines = ["environment\tsoftware\tversion\tbuild\tsource"]
    version_lines.extend(
        (
            f"pipeline\tMETAHICT\t{VERSION}\t-\tmetahict_manager.py",
            f"pipeline\tNextflow\t{bundled_nextflow_version()}\t-\tnextflow/bin/nextflow",
            f"pipeline\tPython\t{platform.python_version()}\t-\truntime",
        )
    )
    for name in ENVIRONMENTS:
        lock_file = LOCK_DIR / f"{name}.explicit.txt"
        prefix = environment_prefix(name)
        lock_lines.append(f"{name}\t{prefix}\t{lock_file}\t{sha256_file(lock_file)}")
        for url in lock_urls(lock_file):
            package = locked_package(url)
            if package is None:
                continue
            package_name, package_version, package_build = package
            version_lines.append(
                f"{name}\t{package_name}\t{package_version}\t{package_build}\texplicit-lock"
            )
    pip_requirements = PROJECT_ROOT / "installation" / "pip-requirements.txt"
    for line in pip_requirements.read_text().splitlines():
        requirement = pinned_python_requirement(line.strip())
        if requirement is None or line.lstrip().startswith("#"):
            continue
        package_name, package_version, package_build = requirement
        version_lines.append(
            f"metahict_env\t{package_name}\t{package_version}\t{package_build}\tpip-requirements"
        )
    (run_dir / "environment_locks.tsv").write_text("\n".join(lock_lines) + "\n")
    (run_dir / "software_versions.tsv").write_text("\n".join(version_lines) + "\n")


def serializable_parameters(args: argparse.Namespace) -> dict[str, object]:
    return {
        key: value
        for key, value in vars(args).items()
        if key != "handler" and isinstance(value, (str, int, float, bool, type(None)))
    }


def initialize_run_metadata(
    args: argparse.Namespace,
    run_id: str,
    run_dir: Path,
    command: Sequence[str],
    previous_run: str | None,
) -> dict[str, object]:
    samplesheet = Path(args.samplesheet).expanduser().resolve()
    config = Path(getattr(args, "config", None) or DEFAULT_CONFIG).expanduser().resolve()
    started_at = datetime.now(timezone.utc)
    parameters = serializable_parameters(args)
    parameters.update(
        {
            "run_id": run_id,
            "resolved_outdir": str(Path(args.outdir).expanduser().resolve()),
            "resolved_report_directory": str(run_dir),
            "resolved_work_directory": str(work_root(args)),
            "nextflow_command": list(command),
        }
    )
    atomic_write_json(run_dir / "parameters.json", parameters)
    sample_digest = sha256_file(samplesheet) if samplesheet.is_file() else None
    digest_record = sample_digest or "MISSING"
    (run_dir / "samplesheet.sha256").write_text(f"{digest_record}  {samplesheet}\n")
    config_digest = sha256_file(config) if config.is_file() else None
    config_record = config_digest or "MISSING"
    (run_dir / "configuration.sha256").write_text(f"{config_record}  {config}\n")
    write_environment_metadata(run_dir)
    provenance: dict[str, object] = {
        "schema_version": 1,
        "run_id": run_id,
        "status": "running",
        "exit_code": None,
        "started_at": started_at.isoformat(),
        "finished_at": None,
        "duration_seconds": None,
        "metahict_version": VERSION,
        "nextflow_version": bundled_nextflow_version(),
        "git_revision": git_revision(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "profile": "local",
        "entry_module": parse_entry_module(getattr(args, "entry_module", "all")),
        "resumed": bool(args.resume),
        "resumes_run_id": previous_run if args.resume else None,
        "samplesheet": str(samplesheet),
        "samplesheet_sha256": sample_digest,
        "configuration": str(config),
        "configuration_sha256": config_digest,
        "outdir": str(Path(args.outdir).expanduser().resolve()),
        "report_directory": str(run_dir),
        "work_directory": str(work_root(args)),
        "command": list(command),
        "command_display": display_command(command),
    }
    atomic_write_json(run_dir / "provenance.json", provenance)
    return provenance


def finalize_provenance(
    path: Path,
    provenance: dict[str, object],
    *,
    status: str,
    exit_code: int,
    started_monotonic: float,
) -> None:
    provenance.update(
        {
            "status": status,
            "exit_code": exit_code,
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "duration_seconds": round(time.monotonic() - started_monotonic, 3),
        }
    )
    atomic_write_json(path, provenance)


def task_work_directory(work: Path, task_hash: str) -> Path | None:
    if "/" not in task_hash:
        return None
    prefix, remainder = task_hash.split("/", 1)
    parent = work / prefix
    if not parent.is_dir():
        return None
    matches = sorted(parent.glob(f"{remainder}*"))
    return matches[-1] if matches else None


def safe_log_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return (cleaned or "process")[:120]


def archive_process_logs(run_dir: Path, work: Path) -> list[dict[str, str]]:
    trace = run_dir / "trace.txt"
    archive = run_dir / "process_logs"
    archive.mkdir(exist_ok=True)
    archived: list[dict[str, str]] = []
    if trace.is_file():
        with trace.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for position, row in enumerate(rows, start=1):
            task_hash = row.get("hash", "")
            reported_work_directory = row.get("workdir", "")
            reported_source = Path(reported_work_directory) if reported_work_directory else None
            source = (
                reported_source
                if reported_source is not None and reported_source.is_dir()
                else task_work_directory(work, task_hash)
            )
            task_name = row.get("name", "process")
            destination = archive / (
                f"{position:04d}-{safe_log_name(task_name)}-{task_hash.replace('/', '-') or 'no-hash'}"
            )
            destination.mkdir(exist_ok=True)
            if source is not None:
                for filename in (
                    ".command.sh",
                    ".command.out",
                    ".command.err",
                    ".command.log",
                    ".exitcode",
                ):
                    candidate = source / filename
                    if candidate.is_file():
                        shutil.copy2(candidate, destination / filename.lstrip("."))
            archived.append(
                {
                    "task_id": row.get("task_id", ""),
                    "name": task_name,
                    "status": row.get("status", ""),
                    "exit": row.get("exit", ""),
                    "hash": task_hash,
                    "work_directory": str(source) if source else "not-found",
                    "log_directory": str(destination),
                }
            )
    index_fields = (
        "task_id",
        "name",
        "status",
        "exit",
        "hash",
        "work_directory",
        "log_directory",
    )
    with (archive / "index.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(archived)
    return archived


def write_failure_summary(
    run_dir: Path, exit_code: int, tasks: Sequence[dict[str, str]]
) -> Path:
    failed = [task for task in tasks if task["status"].upper() in {"FAILED", "ABORTED"}]
    lines = [
        "METAHICT workflow failure summary",
        f"Exit code: {exit_code}",
        f"Run log: {run_dir / 'run.log'}",
        f"Trace: {run_dir / 'trace.txt'}",
    ]
    if failed:
        lines.append("Failed tasks:")
        for task in failed:
            lines.extend(
                (
                    f"  - {task['name']} (status={task['status']}, exit={task['exit']}, hash={task['hash']})",
                    f"    Process logs: {task['log_directory']}",
                    f"    Work directory: {task['work_directory']}",
                )
            )
    else:
        lines.append("No failed task row was available; inspect run.log and nextflow.log.")
    summary = run_dir / "failure_summary.txt"
    summary.write_text("\n".join(lines) + "\n")
    return summary


def preserve_nextflow_engine_log(run_dir: Path) -> None:
    engine_log = PROJECT_ROOT / ".nextflow.log"
    if engine_log.is_file():
        shutil.copy2(engine_log, run_dir / "nextflow.log")


def build_nextflow_command(
    args: argparse.Namespace,
    *,
    report_directory: Path | None = None,
    resource_limits: tuple[int, int] | None = None,
) -> list[str]:
    samplesheet = Path(args.samplesheet).expanduser().resolve()
    config = Path(getattr(args, "config", None) or DEFAULT_CONFIG).expanduser().resolve()
    out_root = Path(args.outdir).expanduser().resolve()
    reports = report_directory or report_root(args)
    work = work_root(args)
    databases = default_database_paths()
    entry_module = parse_entry_module(getattr(args, "entry_module", "all"))
    command = [
        str(PROJECT_ROOT / "nextflow" / "bin" / "nextflow"),
        "run",
        str(PROJECT_ROOT / "nextflow" / "main_dsl2.nf"),
        "-params-file",
        str(config),
        "-profile",
        "local",
        "--entry_module",
        entry_module,
        "--samplesheet",
        str(samplesheet),
        "--out_root",
        str(out_root),
        "--report_dir",
        str(reports),
        "-work-dir",
        str(work),
        "-ansi-log",
        "false",
    ]
    cpu_limit, memory_limit = resource_limits or local_resource_limits()
    memory_limit_mb = max(1, memory_limit // (1024 ** 2))
    command.extend(
        [
            "--local_resource_cpus",
            str(cpu_limit),
            "--local_resource_memory",
            f"{memory_limit_mb} MB",
        ]
    )
    if getattr(args, "threads", None) is not None:
        command.extend(["--threads", str(args.threads)])
    if getattr(args, "memory", None) is not None:
        command.extend(["--memory", args.memory])
    for database_name in required_databases_for_run(args):
        argument_name, nextflow_name, _ = DATABASE_ARGUMENTS[database_name]
        configured = getattr(args, argument_name, None)
        path = Path(configured).expanduser().resolve() if configured else databases[database_name]
        command.extend([f"--{nextflow_name}", str(path)])
    for argument_name, nextflow_name in REUSABLE_STAGE_INPUTS:
        value = getattr(args, argument_name, None)
        if value:
            command.extend(
                [
                    f"--{nextflow_name}",
                    str(Path(value).expanduser().resolve()),
                ]
            )
    mag_dir = getattr(args, "mag_dir", None)
    if mag_dir:
        command.extend(
            [
                "--annotation_mag_dir",
                str(Path(mag_dir).expanduser().resolve()),
            ]
        )
    host_dir = getattr(args, "host_dir", None)
    if host_dir:
        command.extend(
            [
                "--mge_host_dir",
                str(Path(host_dir).expanduser().resolve()),
            ]
        )
    if getattr(args, "scaffolding_skip_checkm2", False):
        command.append("--scaffolding_skip_checkm2")
    if args.resume:
        command.append("-resume")
    return command


def validate_preprocessing_input(
    args: argparse.Namespace,
    samplesheet: Path,
    entry_module: str,
) -> None:
    configured = getattr(args, "preprocessing_dir", None)
    if not configured:
        return
    required_libraries: tuple[str, ...] = {
        "assembly": ("sg",),
        "alignment": ("hic",),
        "coverage": ("sg",),
        "reassembly": ("sg", "hic"),
        "scaffolding": ("hic",),
    }.get(entry_module, ())
    if entry_module == "mge" and not (
        getattr(args, "mge_alignment_dir", None)
        or getattr(args, "mge_contact_dir", None)
    ):
        required_libraries = ("hic",)
    if not required_libraries:
        return

    with samplesheet.open(newline="") as handle:
        rows = [
            row for row in csv.DictReader(handle) if (row.get("sample") or "").strip()
        ]
    sample_libraries = [
        (row["sample"].strip(), library)
        for row in rows
        for library in required_libraries
        if not (library == "sg" and (row.get("long_read_type") or "").strip())
    ]
    if not sample_libraries:
        return
    samples = sorted({sample for sample, _ in sample_libraries})
    configured_text = str(configured)
    if len(samples) > 1 and "{sample}" not in configured_text:
        raise MetahictError(
            "A multi-sample sheet requires {sample} in --preprocessing-dir"
        )
    roots = (
        [
            Path(configured_text.replace("{sample}", sample)).expanduser().resolve()
            for sample in samples
        ]
        if "{sample}" in configured_text
        else [Path(configured_text).expanduser().resolve()]
    )
    missing_roots = [str(root) for root in roots if not root.is_dir()]
    if missing_roots:
        raise MetahictError(
            "Preprocessing output directory not found:\n  "
            + "\n  ".join(missing_roots)
        )
    root_by_sample = (
        {
            sample: Path(configured_text.replace("{sample}", sample))
            .expanduser()
            .resolve()
            for sample in samples
        }
        if "{sample}" in configured_text
        else {sample: roots[0] for sample in samples}
    )
    missing_libraries = [
        str(root_by_sample[sample] / library)
        for sample, library in sample_libraries
        if not (root_by_sample[sample] / library).is_dir()
    ]
    if missing_libraries:
        raise MetahictError(
            "Required preprocessing library directory not found:\n  "
            + "\n  ".join(missing_libraries)
        )


def validate_run_inputs(args: argparse.Namespace) -> None:
    samplesheet = Path(args.samplesheet).expanduser().resolve()
    if not samplesheet.is_file():
        raise MetahictError(f"Samplesheet not found: {samplesheet}")
    config = Path(getattr(args, "config", None) or DEFAULT_CONFIG).expanduser().resolve()
    if not config.is_file():
        raise MetahictError(f"Configuration file not found: {config}")
    validate_configuration_schema(config)
    entry_module = parse_entry_module(getattr(args, "entry_module", "all"))
    with samplesheet.open(newline="") as handle:
        samplesheet_reader = csv.DictReader(handle)
        samplesheet_fields = samplesheet_reader.fieldnames or []
        samplesheet_rows = list(samplesheet_reader)
    if "sample" not in samplesheet_fields:
        raise MetahictError("The samplesheet requires a 'sample' column")
    empty_sample_rows = [
        str(line_number)
        for line_number, row in enumerate(samplesheet_rows, start=2)
        if not (row.get("sample") or "").strip()
    ]
    if empty_sample_rows:
        raise MetahictError(
            "Missing sample identifier for samplesheet row(s) "
            + ", ".join(empty_sample_rows)
        )
    read_entries = {"all", "preprocessing", "assembly", "coverage", "reassembly"}
    if entry_module in read_entries:
        required_columns = {"sg1", "sg2", "hic1", "hic2"}
        missing_columns = sorted(required_columns.difference(samplesheet_fields))
        if missing_columns:
            raise MetahictError(
                "The samplesheet is missing required read column(s): "
                + ", ".join(missing_columns)
            )
        read_errors = []
        for line_number, row in enumerate(samplesheet_rows, start=2):
            long_read_type = (row.get("long_read_type") or "").strip()
            if long_read_type and long_read_type not in LONG_READ_TYPES:
                read_errors.append(
                    f"row {line_number}: invalid long_read_type '{long_read_type}'"
                )
            if not (row.get("sg1") or "").strip():
                read_errors.append(f"row {line_number}: sg1 is required")
            if long_read_type:
                if (row.get("sg2") or "").strip():
                    read_errors.append(
                        f"row {line_number}: sg2 must be empty in long-read mode"
                    )
                if entry_module == "reassembly":
                    read_errors.append(
                        f"row {line_number}: reassembly is not available for long-read samples"
                    )
            elif not (row.get("sg2") or "").strip():
                read_errors.append(
                    f"row {line_number}: sg2 is required unless long_read_type is set"
                )
            for field in ("hic1", "hic2"):
                if not (row.get(field) or "").strip():
                    read_errors.append(f"row {line_number}: {field} is required")
        if read_errors:
            raise MetahictError(
                "Invalid shotgun/Hi-C samplesheet input:\n  "
                + "\n  ".join(read_errors)
            )
    enzyme_required = entry_module in ENZYME_ENTRY_MODULES and not (
        entry_module == "mge" and getattr(args, "mge_contact_dir", None)
    )
    if enzyme_required:
        with samplesheet.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "enzyme" not in reader.fieldnames:
                raise MetahictError(
                    f"The {entry_module} entry requires an 'enzyme' column in the samplesheet"
                )
            missing_enzyme_rows = [
                str(line_number)
                for line_number, row in enumerate(reader, start=2)
                if not (row.get("enzyme") or "").strip()
            ]
        if missing_enzyme_rows:
            raise MetahictError(
                "Missing restriction enzyme for samplesheet row(s) "
                + ", ".join(missing_enzyme_rows)
            )
    with samplesheet.open(newline="") as handle:
        reader = csv.DictReader(handle)
        obsolete_values = [
            f"row {line_number}: {column}"
            for line_number, row in enumerate(reader, start=2)
            for column, value in row.items()
            if column and column.endswith("_extra_args") and (value or "").strip()
        ]
    if obsolete_values:
        raise MetahictError(
            "Samplesheet *_extra_args values are no longer accepted; move these "
            "parameters to the YAML configuration:\n  "
            + "\n  ".join(obsolete_values)
        )
    has_short_read_samples = any(
        not (row.get("long_read_type") or "").strip()
        for row in samplesheet_rows
    )
    required_upstream_options = {
        "assembly": (
            (("preprocessing_dir", "--preprocessing-dir"),)
            if has_short_read_samples
            else ()
        ),
        "alignment": (
            ("assembly_dir", "--assembly-dir"),
            ("preprocessing_dir", "--preprocessing-dir"),
        ),
        "coverage": (
            (("assembly_dir", "--assembly-dir"), ("preprocessing_dir", "--preprocessing-dir"))
            if has_short_read_samples
            else (("assembly_dir", "--assembly-dir"),)
        ),
        "contact": (
            ("assembly_dir", "--assembly-dir"),
            ("alignment_dir", "--alignment-dir"),
        ),
        "binning": (
            ("assembly_dir", "--assembly-dir"),
            ("alignment_dir", "--alignment-dir"),
        ),
        "reassembly": (
            ("binning_dir", "--binning-dir"),
            ("assembly_dir", "--assembly-dir"),
            ("alignment_dir", "--alignment-dir"),
            ("preprocessing_dir", "--preprocessing-dir"),
        ),
        "scaffolding": (("preprocessing_dir", "--preprocessing-dir"),),
    }
    missing_upstream_options = [
        option
        for attribute, option in required_upstream_options.get(entry_module, ())
        if not getattr(args, attribute, None)
    ]
    if missing_upstream_options:
        raise MetahictError(
            f"The {entry_module} entry requires explicit upstream input(s): "
            + ", ".join(missing_upstream_options)
        )
    if entry_module == "annotation":
        configured_mag_dir = getattr(args, "mag_dir", None)
        if not configured_mag_dir:
            raise MetahictError(
                "The annotation entry requires a directory containing MAG FASTAs: "
                "--mag-dir PATH"
            )
        with samplesheet.open(newline="") as handle:
            samples = [
                row.get("sample", "").strip()
                for row in csv.DictReader(handle)
                if row.get("sample", "").strip()
            ]
        if not samples:
            raise MetahictError(
                "The annotation samplesheet requires at least one non-empty sample row"
            )
        configured_text = str(configured_mag_dir)
        if len(samples) > 1 and "{sample}" not in configured_text:
            raise MetahictError(
                "A multi-sample annotation sheet requires {sample} in --mag-dir "
                "so each sample resolves to a separate MAG directory"
            )
        if "{sample}" in configured_text:
            mag_directories = [
                Path(configured_text.replace("{sample}", sample)).expanduser().resolve()
                for sample in samples
            ]
        else:
            mag_directories = [Path(configured_text).expanduser().resolve()]
        missing_directories = [
            str(directory) for directory in mag_directories if not directory.is_dir()
        ]
        if missing_directories:
            raise MetahictError(
                "Annotation MAG directory not found:\n  "
                + "\n  ".join(missing_directories)
            )
        supported_suffixes = {".fa", ".fasta", ".fna"}
        empty_directories = [
            str(directory)
            for directory in mag_directories
            if not any(
                path.is_file() and path.suffix.lower() in supported_suffixes
                for path in directory.iterdir()
            )
        ]
        if empty_directories:
            raise MetahictError(
                "Annotation MAG directory contains no .fa, .fasta, or .fna genomes:\n  "
                + "\n  ".join(empty_directories)
            )
    elif getattr(args, "mag_dir", None):
        raise MetahictError(
            "--mag-dir is available only with --entry-module annotation"
        )
    if entry_module == "scaffolding":
        configured_bin = getattr(args, "scaffolding_bin", None)
        if not configured_bin:
            raise MetahictError(
                "The scaffolding entry requires one explicit bin FASTA: "
                "--scaffolding-bin PATH"
            )
        if getattr(args, "alignment_dir", None):
            raise MetahictError(
                "The scaffolding entry does not use --alignment-dir. "
                "Pass --scaffolding-bin and --preprocessing-dir instead."
            )
        with samplesheet.open(newline="") as handle:
            samples = [
                row.get("sample", "").strip()
                for row in csv.DictReader(handle)
                if row.get("sample", "").strip()
            ]
        configured_text = str(configured_bin)
        if "{sample}" in configured_text:
            bin_paths = [
                Path(configured_text.replace("{sample}", sample)).expanduser().resolve()
                for sample in samples
            ]
        else:
            bin_paths = [Path(configured_text).expanduser().resolve()]
        missing_bins = [str(path) for path in bin_paths if not path.is_file()]
        if missing_bins:
            raise MetahictError(
                "Scaffolding bin FASTA not found:\n  " + "\n  ".join(missing_bins)
            )
        configured_bam = getattr(args, "scaffolding_bam", None)
        if configured_bam:
            configured_bam_text = str(configured_bam)
            if "{sample}" in configured_bam_text:
                bam_paths = [
                    Path(configured_bam_text.replace("{sample}", sample))
                    .expanduser()
                    .resolve()
                    for sample in samples
                ]
            else:
                bam_paths = [Path(configured_bam_text).expanduser().resolve()]
            missing_bams = [str(path) for path in bam_paths if not path.is_file()]
            if missing_bams:
                raise MetahictError(
                    "Scaffolding BAM not found:\n  " + "\n  ".join(missing_bams)
                )
    elif getattr(args, "scaffolding_bam", None):
        raise MetahictError(
            "--scaffolding-bam is available only with --entry-module scaffolding"
        )
    if entry_module == "mge":
        configured_fasta = getattr(args, "mge_fasta", None)
        if not configured_fasta:
            raise MetahictError("The MGE entry requires --fasta PATH")
        with samplesheet.open(newline="") as handle:
            samples = [
                row.get("sample", "").strip()
                for row in csv.DictReader(handle)
                if row.get("sample", "").strip()
            ]
        if not samples:
            raise MetahictError(
                "The MGE samplesheet requires at least one non-empty sample row"
            )
        fasta_text = str(configured_fasta)
        if len(samples) > 1 and "{sample}" not in fasta_text:
            raise MetahictError(
                "A multi-sample MGE sheet requires {sample} in --fasta"
            )
        fasta_paths = (
            [
                Path(fasta_text.replace("{sample}", sample)).expanduser().resolve()
                for sample in samples
            ]
            if "{sample}" in fasta_text
            else [Path(fasta_text).expanduser().resolve()]
        )
        missing_fastas = [str(path) for path in fasta_paths if not path.is_file()]
        if missing_fastas:
            raise MetahictError(
                "MGE FASTA not found:\n  " + "\n  ".join(missing_fastas)
            )
        configured_host_dir = getattr(args, "host_dir", None)
        if not configured_host_dir:
            raise MetahictError("The MGE entry requires --host-dir PATH")
        host_text = str(configured_host_dir)
        if len(samples) > 1 and "{sample}" not in host_text:
            raise MetahictError(
                "A multi-sample MGE sheet requires {sample} in --host-dir"
            )
        host_directories = (
            [
                Path(host_text.replace("{sample}", sample)).expanduser().resolve()
                for sample in samples
            ]
            if "{sample}" in host_text
            else [Path(host_text).expanduser().resolve()]
        )
        missing_host_directories = [
            str(directory) for directory in host_directories if not directory.is_dir()
        ]
        if missing_host_directories:
            raise MetahictError(
                "MGE host directory not found:\n  "
                + "\n  ".join(missing_host_directories)
            )
        supported_suffixes = {".fa", ".fasta", ".fna"}
        empty_host_directories = [
            str(directory)
            for directory in host_directories
            if not any(
                path.is_file() and path.suffix.lower() in supported_suffixes
                for path in directory.iterdir()
            )
        ]
        if empty_host_directories:
            raise MetahictError(
                "MGE host directory contains no .fa, .fasta, or .fna genomes:\n  "
                + "\n  ".join(empty_host_directories)
            )
        if not getattr(args, "mge_contact_dir", None) and not getattr(
            args, "mge_alignment_dir", None
        ) and not (
            getattr(args, "preprocessing_dir", None)
        ):
            raise MetahictError(
                "The MGE entry requires --preprocessing-dir unless "
                "--mge-alignment-dir or --mge-contact-dir is supplied"
            )
    elif any(
        getattr(args, name, None)
        for name in ("mge_fasta", "host_dir", "mge_alignment_dir", "mge_contact_dir")
    ):
        raise MetahictError(
            "--fasta, --host-dir, --mge-alignment-dir, and --mge-contact-dir are "
            "available only with --entry-module mge"
        )
    validate_preprocessing_input(args, samplesheet, entry_module)
    verify_runtime(verbose=getattr(args, "verbose_preflight", False))
    if getattr(args, "check_outputs", False) and entry_module != "all":
        raise MetahictError(
            "--check-outputs is available only for the complete core example-dataset workflow"
        )
    required_databases = required_databases_for_run(args)
    if not required_databases:
        return
    db_paths = default_database_paths()
    selected = {}
    for database_name in required_databases:
        argument_name, _, label = DATABASE_ARGUMENTS[database_name]
        configured = getattr(args, argument_name, None)
        selected[label] = (
            Path(configured).expanduser().resolve() if configured else db_paths[database_name]
        )
    missing = [f"{label}: {path}" for label, path in selected.items() if not path.exists()]
    if missing:
        raise MetahictError("Missing database path(s):\n  " + "\n  ".join(missing))


def command_run(args: argparse.Namespace) -> None:
    resource_limits = local_resource_limits()
    if args.show_command:
        command = build_nextflow_command(args, resource_limits=resource_limits)
        warn_if_resources_are_capped(args, *resource_limits)
        print("[INFO] Nextflow command:")
        print(display_command(command))
        return

    reports = report_root(args)
    reports.mkdir(parents=True, exist_ok=True)
    run_id, run_dir = create_run_directory(reports)
    previous_run = prepare_report_links(reports, run_id, run_dir)
    command = build_nextflow_command(
        args, report_directory=run_dir, resource_limits=resource_limits
    )
    provenance = initialize_run_metadata(args, run_id, run_dir, command, previous_run)
    log_path = run_dir / "run.log"
    started_monotonic = time.monotonic()
    terminal_stdout = sys.stdout
    terminal_stderr = sys.stderr
    status = 1
    failure: BaseException | None = None

    with log_path.open("x", buffering=1) as log_handle:
        with contextlib.redirect_stdout(TeeTextIO(terminal_stdout, log_handle)), contextlib.redirect_stderr(
            TeeTextIO(terminal_stderr, log_handle)
        ):
            print(f"[INFO] METAHICT run ID: {run_id}")
            print(f"[INFO] Run directory: {run_dir}")
            print("[INFO] Nextflow command:")
            print(display_command(command))
            try:
                validate_run_inputs(args)
                warn_if_resources_are_capped(args, *resource_limits)
                process = subprocess.Popen(
                    command,
                    cwd=str(PROJECT_ROOT),
                    env=nextflow_environment(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                assert process.stdout is not None
                for line in process.stdout:
                    sys.stdout.write(line)
                status = process.wait()
                if status:
                    raise MetahictError(f"Nextflow failed with exit code {status}")
                if args.check_outputs:
                    check_command = [
                        sys.executable,
                        PROJECT_ROOT / "nextflow" / "bin" / "check_expected_outputs.py",
                        "--root",
                        Path(args.outdir).expanduser().resolve(),
                        "--manifest",
                        PROJECT_ROOT / "nextflow" / "tests" / "expected" / "example_dataset_outputs.tsv",
                    ]
                    print(f"[RUN] {display_command(check_command)}", flush=True)
                    checker = subprocess.Popen(
                        [str(item) for item in check_command],
                        cwd=str(PROJECT_ROOT),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                    )
                    assert checker.stdout is not None
                    for line in checker.stdout:
                        sys.stdout.write(line)
                    check_status = checker.wait()
                    if check_status:
                        status = check_status
                        raise MetahictError(
                            f"Expected-output validation failed with exit code {check_status}"
                        )
                status = 0
            except BaseException as error:
                failure = error
                if isinstance(error, subprocess.CalledProcessError):
                    status = error.returncode or 1
                elif isinstance(error, KeyboardInterrupt):
                    status = 130
                elif status == 0:
                    status = 1
                log_handle.write(f"[ERROR] {error}\n")
            finally:
                try:
                    preserve_nextflow_engine_log(run_dir)
                    archived_tasks = archive_process_logs(run_dir, work_root(args))
                except Exception as error:
                    archived_tasks = []
                    print(f"[WARN] Could not archive all process logs: {error}", file=sys.stderr)
                if status:
                    summary = write_failure_summary(run_dir, status, archived_tasks)
                    print(f"[INFO] Failure summary: {summary}")
                finalize_provenance(
                    run_dir / "provenance.json",
                    provenance,
                    status="completed" if status == 0 else "failed",
                    exit_code=status,
                    started_monotonic=started_monotonic,
                )
                print(f"[INFO] Run log: {log_path}")
                print(f"[INFO] Process-log index: {run_dir / 'process_logs' / 'index.tsv'}")
                if status == 0:
                    print(f"[PASS] METAHICT workflow completed; immutable report directory: {run_dir}")

    if failure is not None:
        raise failure


def command_config(args: argparse.Namespace) -> None:
    output = Path(args.output).expanduser().resolve()
    if output.exists() and not args.force:
        raise MetahictError(
            f"Configuration already exists: {output}. Use --force to replace it."
        )
    text = DEFAULT_CONFIG.read_text()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text)
    print(f"[PASS] Configuration created: {output}")


def command_samplesheet(args: argparse.Namespace) -> None:
    template = PROJECT_ROOT / "nextflow" / "assets" / "samplesheet_template.csv"
    with template.open(newline="") as handle:
        fields = next(csv.reader(handle))
    long_read_type = (getattr(args, "long_read_type", None) or "").strip()
    sg_r2 = getattr(args, "sg_r2", None)
    if long_read_type and sg_r2:
        raise MetahictError(
            "--sg-r2 must be omitted when --long-read-type is supplied"
        )
    if not long_read_type and not sg_r2:
        raise MetahictError(
            "--sg-r2 is required for paired short reads; alternatively supply "
            "--long-read-type for a single long-read shotgun file"
        )
    paths = {
        "sg1": Path(args.sg_r1).expanduser().resolve(),
        "sg2": Path(sg_r2).expanduser().resolve() if sg_r2 else None,
        "hic1": Path(args.hic_r1).expanduser().resolve(),
        "hic2": Path(args.hic_r2).expanduser().resolve(),
    }
    if not args.allow_missing:
        missing = [
            str(path) for path in paths.values() if path is not None and not path.is_file()
        ]
        if missing:
            raise MetahictError(
                "Input read file(s) not found:\n  " + "\n  ".join(missing)
            )
    row = {field: "" for field in fields}
    row.update(
        {
            "sample": args.sample,
            "sg1": str(paths["sg1"]),
            "sg2": str(paths["sg2"]) if paths["sg2"] else "",
            "hic1": str(paths["hic1"]),
            "hic2": str(paths["hic2"]),
            "enzyme": args.enzyme,
            "long_read_type": long_read_type,
        }
    )
    output = Path(args.output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    append_existing = args.append and output.is_file()
    if append_existing:
        with output.open(newline="") as handle:
            existing_fields = next(csv.reader(handle), [])
        legacy_fields = [field for field in fields if field != "long_read_type"]
        if existing_fields == legacy_fields and not long_read_type:
            fields = legacy_fields
            row.pop("long_read_type", None)
        elif existing_fields != fields:
            raise MetahictError(f"Existing samplesheet has an incompatible header: {output}")
    with output.open("a" if append_existing else "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if not append_existing:
            writer.writeheader()
        writer.writerow(row)
    print(f"[PASS] Samplesheet created: {output}")


def _plain_help_line(value: str) -> str:
    """Remove lightweight Markdown markup for terminal help output."""
    value = re.sub(r"\[([^]]+)\]\([^)]+\)", r"\1", value)
    return value.replace("`", "").replace("**", "").strip()


def configured_module_keys(module: str) -> set[str]:
    """Return scalar key paths below one `modules.<name>` template block."""
    keys: set[str] = set()
    in_modules = False
    in_selected = False
    nested: list[tuple[int, str]] = []
    for raw_line in DEFAULT_CONFIG.read_text().splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip())
        if indent == 0:
            in_modules = stripped == "modules:"
            in_selected = False
            continue
        if not in_modules:
            continue
        if indent == 2 and stripped.endswith(":"):
            in_selected = stripped[:-1] == module
            nested = []
            continue
        if not in_selected or indent <= 2 or ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        while nested and indent <= nested[-1][0]:
            nested.pop()
        path = [item[1] for item in nested] + [key]
        if value.strip():
            keys.add(".".join(path))
        else:
            nested.append((indent, key))
    return keys


def _markdown_help_to_text(markdown: str, module: str) -> str:
    """Render the maintained module manual as readable terminal text."""
    output: list[str] = []
    in_code = False
    table_header: list[str] | None = None
    current_section = ""
    module_keys = configured_module_keys(module)
    for raw_line in markdown.splitlines():
        line = raw_line.rstrip()
        if line.strip().startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            output.append(f"  {line}" if line else "")
            continue
        if re.match(r"^\|(?:\s*:?-+:?\s*\|)+$", line):
            continue
        if line.startswith("|") and line.endswith("|"):
            cells = [_plain_help_line(cell) for cell in line.strip("|").split("|")]
            if cells and cells[0] in {"Parameter", "Output"}:
                table_header = cells
                continue
            if len(cells) >= 3:
                candidate = cells[0].removeprefix("--").replace("-", "_")
                if candidate in module_keys:
                    cells[0] = f"modules.{module}.{candidate}"
                elif candidate.startswith("no_") and candidate[3:] in module_keys:
                    cells[0] = f"modules.{module}.{candidate[3:]} = false"
                output.append(f"  {cells[0]}")
                if table_header and len(table_header) > 1 and cells[1]:
                    output.append(f"      {table_header[1]}: {cells[1]}")
                output.append(f"      {cells[2]}")
            elif len(cells) == 2:
                output.append(f"  {cells[0]}: {cells[1]}")
            continue
        table_header = None
        heading = re.match(r"^(#{1,6})\s+(.+)$", line)
        if heading:
            title = _plain_help_line(heading.group(2))
            current_section = title.lower()
            output.extend(["", title.upper(), ""])
            continue
        output.append(_plain_help_line(line) if line else "")
    return "\n".join(output).strip()


def default_module_resource(module: str) -> tuple[str, str]:
    """Read the distributed resource row without requiring a YAML package."""
    pattern = re.compile(
        rf"^\s{{2}}{re.escape(module)}:\s*\{{threads:\s*([^,]+),\s*memory:\s*[\"']([^\"']+)[\"']\}}\s*$",
        re.MULTILINE,
    )
    match = pattern.search(DEFAULT_CONFIG.read_text())
    if match is None:
        raise MetahictError(f"Missing distributed resource defaults for module: {module}")
    return match.group(1).strip(), match.group(2).strip()


def module_help_text(module: str) -> str:
    """Return detailed help for one public Nextflow entry module."""
    document = MODULE_HELP_DOCS[module]
    if not document.is_file():
        raise MetahictError(f"Module help document is missing: {document}")
    threads, memory = default_module_resource(module)
    resource_behavior = MODULE_RESOURCE_BEHAVIOR[module]
    yaml_sections = f"modules.{module}"
    preamble = f"""usage: ./metahict run --entry-module {module} --samplesheet FILE --config FILE --outdir DIR [options]

MODULE-SPECIFIC HELP

The options that select files or databases are ./metahict run command-line
arguments. Algorithm settings are YAML keys under {yaml_sections}; do not put
those YAML keys directly on the command line.

RESOURCE BEHAVIOR

Default: resources.{module}.threads={threads}, resources.{module}.memory={memory}
{resource_behavior}

-t/--threads and -m/--memory override these configured values for every selected
module in one run. When neither option is supplied, the YAML values are used.
If a request exceeds detected local capacity, Nextflow warns and caps it. The
effective values are passed to tools and recorded in the Nextflow trace.
"""
    return f"{preamble.rstrip()}\n\n{_markdown_help_to_text(document.read_text(), module)}\n"


def requested_module_help(argv: Sequence[str]) -> str | None:
    """Recognize `run --entry-module NAME -h/--help` before required parsing."""
    arguments = list(argv)
    if not arguments or arguments[0] != "run" or not ({"-h", "--help"} & set(arguments)):
        return None
    selected: str | None = None
    for index, value in enumerate(arguments):
        if value == "--entry-module" and index + 1 < len(arguments):
            selected = arguments[index + 1]
        elif value.startswith("--entry-module="):
            selected = value.split("=", 1)[1]
    if selected is None:
        return None
    normalized = selected.strip().lower()
    return normalized if normalized in MODULE_HELP_DOCS else None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="./metahict",
        description=(
            "METAHICT is a Nextflow DSL2 workflow for genome-resolved microbiome "
            "analysis with short- or long-read shotgun metagenomic data and "
            "paired-end metagenomic Hi-C reads.\n\n"
            "Use this command to install the locked runtime and databases, create "
            "the configuration and samplesheet, validate the installation, and run "
            "the complete workflow or one selected module."
        ),
        epilog=(
            "GETTING STARTED\n\n"
            "  1. Check the host and install the reproducible Linux runtime:\n"
            "     ./metahict doctor\n"
            "     ./metahict install\n\n"
            "  2. Test the core workflow and standalone scaffolding entry:\n"
            "     ./metahict test workflow\n\n"
            "  3. Install the required reference databases:\n"
            "     ./metahict database all\n\n"
            "  4. Check the installed runtime and databases:\n"
            "     ./metahict doctor --runtime --databases\n\n"
            "  5. Create the editable configuration and inspect samplesheet inputs:\n"
            "     ./metahict config\n"
            "     ./metahict samplesheet --help\n\n"
            "  6. Inspect run options, then run the workflow:\n"
            "     ./metahict run --help\n"
            "     ./metahict run --samplesheet samplesheet.csv \\\n"
            "       --config metahict_configuration.yaml --outdir results\n\n"
            "DETAILED HELP\n\n"
            "  ./metahict COMMAND --help\n"
            "  ./metahict run --entry-module MODULE --help\n\n"
            "Supported modules: preprocessing, assembly, alignment, coverage,\n"
            "contact, binning, reassembly, scaffolding, annotation, and mge."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"METAHICT {VERSION}")
    subparsers = parser.add_subparsers(
        dest="command", required=True, title="commands", metavar="COMMAND"
    )

    install_parser = subparsers.add_parser(
        "install", help="Install the exact Linux dependency bundle"
    )
    install_parser.set_defaults(handler=command_install)

    database_parser = subparsers.add_parser(
        "database", help="Download one or all required reference databases"
    )
    database_parser.add_argument(
        "target", choices=("all", "checkm", "checkm2", "gtdbtk", "genomad")
    )
    database_parser.add_argument("--path", help="Destination for one selected database")
    database_parser.add_argument("--root", help="Parent directory used with target 'all'")
    database_parser.set_defaults(handler=command_database)

    doctor_parser = subparsers.add_parser(
        "doctor", help="Check the host, locks, environments, and databases"
    )
    doctor_parser.add_argument("--runtime", action="store_true", help="Verify all installed environments")
    doctor_parser.add_argument("--databases", action="store_true", help="Verify default database paths")
    doctor_parser.set_defaults(handler=command_doctor)

    verify_parser = subparsers.add_parser(
        "verify", help="Verify the exact installed runtime used by Nextflow"
    )
    verify_parser.set_defaults(handler=lambda _: verify_runtime())

    test_parser = subparsers.add_parser(
        "test",
        help="Run source, workflow, or bundled-example tests",
        description=(
            "Run developer source checks, fast stub workflow checks, or the "
            "bundled real-data example. The workflow test exercises the default "
            "complete workflow and the standalone scaffolding entry separately."
        ),
    )
    test_parser.add_argument(
        "scope",
        nargs="?",
        choices=("source", "workflow", "example", "all"),
        default="all",
    )
    test_parser.add_argument(
        "--outdir",
        default="results",
        help="Output directory used by `test example`",
    )
    test_parser.add_argument(
        "--config",
        default=str(
            PROJECT_ROOT
            / "nextflow"
            / "assets"
            / "example_dataset_configuration.yaml"
        ),
        help="Configuration used by `test example`",
    )
    test_parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume cached tasks during `test example`",
    )
    test_parser.add_argument(
        "--checkm-db", help="CheckM reference directory for `test example`"
    )
    test_parser.add_argument(
        "--checkm2-db", help="CheckM2 database for `test example`"
    )
    test_parser.add_argument(
        "--gtdbtk-db", help="GTDB-Tk reference directory for `test example`"
    )
    test_parser.add_argument(
        "--genomad-db", help="geNomad reference directory for `test example`"
    )
    test_parser.set_defaults(handler=command_test)

    compare_parser = subparsers.add_parser(
        "compare", help="Compare candidate scientific results with an accepted baseline"
    )
    compare_parser.add_argument("--baseline", required=True, help="Accepted sample-result directory")
    compare_parser.add_argument("--candidate", required=True, help="Candidate sample-result directory")
    compare_parser.add_argument("--manifest", help="Comparison manifest; defaults to the release manifest")
    compare_parser.add_argument("--report", default="validation/scientific-regression.tsv")
    compare_parser.set_defaults(handler=command_compare)

    config_parser = subparsers.add_parser(
        "config", help="Create an editable all-module YAML configuration"
    )
    config_parser.add_argument("--output", default="metahict_configuration.yaml")
    config_parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    config_parser.set_defaults(handler=command_config)

    samplesheet_parser = subparsers.add_parser(
        "samplesheet", help="Create a one-sample input and enzyme CSV"
    )
    samplesheet_parser.add_argument("--sample", required=True, help="Short sample identifier")
    samplesheet_parser.add_argument(
        "--sg-r1",
        required=True,
        help=(
            "Shotgun forward FASTQ for paired short reads, or the single "
            "shotgun FASTA/FASTQ file for long reads"
        ),
    )
    samplesheet_parser.add_argument(
        "--sg-r2",
        help="Shotgun reverse FASTQ; required for short reads and omitted for long reads",
    )
    samplesheet_parser.add_argument(
        "--long-read-type",
        choices=LONG_READ_TYPES,
        help="metaFlye input type; setting this selects single-file long-read mode",
    )
    samplesheet_parser.add_argument("--hic-r1", required=True, help="Hi-C forward FASTQ")
    samplesheet_parser.add_argument("--hic-r2", required=True, help="Hi-C reverse FASTQ")
    samplesheet_parser.add_argument(
        "--enzyme",
        required=True,
        help="Hi-C restriction enzyme(s), e.g. Sau3AI,MluCI",
    )
    samplesheet_parser.add_argument("--output", default="samplesheet.csv")
    samplesheet_parser.add_argument(
        "--append", action="store_true", help="Append this sample to an existing generated CSV"
    )
    samplesheet_parser.add_argument(
        "--allow-missing", action="store_true", help="Write paths without checking files"
    )
    samplesheet_parser.set_defaults(handler=command_samplesheet)

    run_parser = subparsers.add_parser(
        "run",
        help="Run the complete workflow or one named Nextflow stage",
        description=(
            "Run the complete METAHICT core workflow (default) or one selected "
            "module. Scaffolding is an optional selected-module analysis and is "
            "not run automatically.\n\n"
            "Read-file paths and restriction enzymes come from --samplesheet. "
            "Algorithm settings and normal per-module resources come from --config. "
            "Explicit command-line -t/--threads and -m/--memory values override "
            "the YAML resources for one run."
        ),
        epilog=(
            "Examples:\n"
            "  # Complete workflow\n"
            "  ./metahict run --samplesheet samplesheet.csv --config metahict_configuration.yaml \\\n"
            "    --outdir results\n\n"
            "  # Detailed inputs, YAML keys, outputs, and command for one module\n"
            "  ./metahict run --entry-module binning --help\n\n"
            "  # Resume an interrupted run with the same inputs and output directory\n"
            "  ./metahict run --samplesheet samplesheet.csv --config metahict_configuration.yaml \\\n"
            "    --outdir results --resume\n\n"
            "Module-specific help:\n"
            "  preprocessing, assembly, alignment, coverage, contact, binning,\n"
            "  reassembly, scaffolding, annotation, and mge"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    core_group = run_parser.add_argument_group("required and workflow selection")
    core_group.add_argument(
        "--samplesheet",
        required=True,
        metavar="FILE",
        help=(
            "CSV with sample, sg1, sg2, hic1, hic2, enzyme, and optional "
            "long_read_type columns"
        ),
    )
    core_group.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG),
        metavar="FILE",
        help="All-module YAML parameters and resources (default: distributed template)",
    )
    core_group.add_argument(
        "--outdir", required=True, metavar="DIR", help="Top-level publication directory"
    )
    core_group.add_argument(
        "--entry-module",
        type=parse_entry_module,
        default="all",
        metavar="NAME",
        help=(
            "Workflow entry: all, preprocessing, assembly, alignment, coverage, "
            "contact, binning, reassembly, scaffolding, annotation, or mge "
            "(default: all)"
        ),
    )

    resource_group = run_parser.add_argument_group("run-wide resource overrides")
    resource_group.add_argument(
        "-t",
        "--threads",
        type=positive_integer,
        default=None,
        help=(
            "Override resources.<stage>.threads for every selected module in "
            "this run; local capacity caps the effective value passed to tools"
        ),
    )
    resource_group.add_argument(
        "-m",
        "--memory",
        type=parse_memory_size,
        default=None,
        metavar="SIZE",
        help=(
            "Override resources.<stage>.memory for every selected module in "
            "this run, e.g. 96G; local capacity caps the effective allocation "
            "used to derive tool-memory limits"
        ),
    )
    database_group = run_parser.add_argument_group("reference-database overrides")
    database_group.add_argument("--checkm-db", metavar="DIR", help="CheckM reference directory")
    database_group.add_argument("--checkm2-db", metavar="FILE", help="CheckM2 DIAMOND database")
    database_group.add_argument("--gtdbtk-db", metavar="DIR", help="GTDB-Tk reference directory")
    database_group.add_argument("--genomad-db", metavar="DIR", help="geNomad reference directory")

    upstream_group = run_parser.add_argument_group("reusable upstream results for selected-module runs")
    upstream_group.add_argument(
        "--preprocessing-dir",
        metavar="DIR",
        help=(
            "Preprocessing output directory containing sg/ and hic/; use "
            "{sample} in the path for a multi-sample sheet"
        ),
    )
    upstream_group.add_argument(
        "--assembly-dir", metavar="DIR", help="Assembly module output directory"
    )
    upstream_group.add_argument(
        "--alignment-dir", metavar="DIR", help="Alignment module output directory"
    )
    upstream_group.add_argument(
        "--binning-dir", metavar="DIR", help="Binning module output directory"
    )
    specialized_group = run_parser.add_argument_group("annotation, scaffolding, and MGE inputs")
    specialized_group.add_argument(
        "--mag-dir",
        dest="mag_dir",
        metavar="DIR",
        help=(
            "Directory of MAG FASTAs required by annotation. "
            "Use {sample} in the path for a multi-sample sheet"
        ),
    )
    specialized_group.add_argument(
        "--scaffolding-bin",
        metavar="FASTA",
        help=(
            "One bin FASTA to scaffold; required for --entry-module scaffolding. "
            "Use {sample} in the path for a multi-sample sheet"
        ),
    )
    specialized_group.add_argument(
        "--scaffolding-bam",
        metavar="BAM",
        help=(
            "Optional BAM aligned to --scaffolding-bin. When omitted, METAHICT "
            "aligns the cleaned Hi-C reads to the selected bin"
        ),
    )
    specialized_group.add_argument(
        "--scaffolding-skip-checkm2",
        action="store_true",
        help="Skip scaffolding quality evaluation and its CheckM2 database requirement",
    )
    specialized_group.add_argument(
        "--fasta",
        dest="mge_fasta",
        metavar="FASTA",
        help=(
            "Metagenome FASTA required by --entry-module mge. Use {sample} in "
            "the path for a multi-sample sheet"
        ),
    )
    specialized_group.add_argument(
        "--host-dir",
        dest="host_dir",
        metavar="DIR",
        help=(
            "Directory of host-genome FASTAs required by --entry-module mge. "
            "Use {sample} in the path for a multi-sample sheet"
        ),
    )
    specialized_group.add_argument(
        "--mge-alignment-dir", metavar="DIR", help="Reuse an MGE-specific alignment directory"
    )
    specialized_group.add_argument(
        "--mge-contact-dir", metavar="DIR", help="Reuse an MGE-specific contact directory"
    )

    control_group = run_parser.add_argument_group("run control, logging, and diagnostics")
    control_group.add_argument(
        "--report-dir", metavar="DIR", help="Override the Nextflow report directory"
    )
    control_group.add_argument(
        "--work-dir", metavar="DIR", help="Override the Nextflow work directory"
    )
    control_group.add_argument(
        "--resume", action="store_true", help="Reuse cached successful tasks from the work directory"
    )
    control_group.add_argument(
        "--verbose-preflight",
        action="store_true",
        help="Print every lock and environment verification instead of one summary line",
    )
    control_group.add_argument(
        "--check-outputs",
        action="store_true",
        help="Check the documented core example-dataset outputs after completion",
    )
    control_group.add_argument(
        "--show-command", action="store_true", help="Print the generated Nextflow command and exit"
    )
    run_parser.set_defaults(handler=command_run)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    selected_help = requested_module_help(arguments)
    if selected_help is not None:
        try:
            print(module_help_text(selected_help), end="")
        except MetahictError as error:
            print(f"[ERROR] {error}", file=sys.stderr)
            return 1
        return 0
    parser = build_parser()
    args = parser.parse_args(arguments)
    try:
        args.handler(args)
    except MetahictError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as error:
        print(f"[ERROR] Command failed with exit code {error.returncode}", file=sys.stderr)
        return error.returncode or 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
