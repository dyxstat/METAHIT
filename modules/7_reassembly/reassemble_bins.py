#!/usr/bin/env python3
"""Recruit reads and reassemble METAHICT bins without shell orchestration.

The bin-selection workflow was adapted from metaWRAP's GPL-3.0 reassembly
module; provenance is recorded in docs/third_party.md. External bioinformatics
programs are installed dependencies and are invoked with explicit argv lists.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Iterator, Optional, Sequence


class ReassemblyError(RuntimeError):
    """Expected reassembly failure."""


SCRIPT_DIR = Path(__file__).resolve().parent
HELPERS = SCRIPT_DIR


def run(
    command: Sequence[object],
    *,
    env: Optional[dict[str, str]] = None,
    capture: bool = False,
    stdout=None,
) -> subprocess.CompletedProcess:
    argv = [str(item) for item in command]
    print(f"[RUN] {shlex.join(argv)}", flush=True)
    return subprocess.run(
        argv,
        check=True,
        env=env,
        text=stdout is None,
        capture_output=capture,
        stdout=stdout,
    )


def start_process(command: Sequence[object], **kwargs) -> subprocess.Popen:
    """Start a pipeline command after normalizing every argv item to text."""
    return subprocess.Popen([str(item) for item in command], **kwargs)


def replace_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def required_program(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise ReassemblyError(f"Required executable not found: {name}")
    return path


def read_fasta(path: Path) -> Iterator[tuple[str, list[str]]]:
    header: Optional[str] = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if stripped.startswith(">"):
                if header is not None:
                    yield header, chunks
                header = stripped[1:].split()[0]
                chunks = []
            elif header is not None:
                chunks.append(stripped)
    if header is not None:
        yield header, chunks


def write_minimum_length_fasta(source: Path, destination: Path, minimum: int) -> int:
    count = 0
    with destination.open("w") as output:
        for identifier, chunks in read_fasta(source):
            sequence = "".join(chunks)
            if len(sequence) < minimum:
                continue
            output.write(f">{identifier}\n{sequence}\n")
            count += 1
    return count


def concatenate_fastas(directory: Path, output: Path) -> int:
    count = 0
    with output.open("wb") as destination:
        for fasta in sorted(path for path in directory.iterdir() if path.is_file()):
            with fasta.open("rb") as source:
                shutil.copyfileobj(source, destination)
            count += 1
    return count


def recruit_reads(
    assembly: Path,
    reads_1: Path,
    reads_2: Path,
    original_bins: Path,
    output: Path,
    threads: int,
    strict_cutoff: int,
    permissive_cutoff: int,
) -> None:
    """Stream one BWA SAM output to read recruitment and unmapped extraction."""
    bwa = required_program("bwa")
    samtools = required_program("samtools")
    if not Path(f"{assembly}.amb").is_file():
        run([bwa, "index", assembly])
    else:
        print(f"[INFO] Reusing existing BWA index: {assembly}")
    reads_output = output / "reads_for_reassembly"
    replace_path(reads_output)
    reads_output.mkdir()
    filter_command = [str(item) for item in [
        sys.executable,
        HELPERS / "filter_reads_for_bin_reassembly.py",
        original_bins,
        reads_output,
        strict_cutoff,
        permissive_cutoff,
    ]]
    print(f"[RUN] streaming {bwa} mem output to samtools and {shlex.join(filter_command)}")
    bwa_process = start_process(
        [bwa, "mem", "-t", str(threads), str(assembly), str(reads_1), str(reads_2)],
        stdout=subprocess.PIPE,
    )
    view_process = start_process(
        [samtools, "view", "-b", "-f", "12", "-"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )
    fastq_process = start_process(
        [
            samtools, "fastq", "-1", str(output / "unmapped_shotgun_1.fastq"),
            "-2", str(output / "unmapped_shotgun_2.fastq"), "-0", os.devnull,
            "-s", os.devnull, "-",
        ],
        stdin=view_process.stdout,
    )
    assert view_process.stdout is not None
    view_process.stdout.close()
    filter_process = start_process(filter_command, stdin=subprocess.PIPE)
    assert bwa_process.stdout is not None
    assert view_process.stdin is not None
    assert filter_process.stdin is not None
    try:
        while True:
            chunk = bwa_process.stdout.read(1024 * 1024)
            if not chunk:
                break
            view_process.stdin.write(chunk)
            filter_process.stdin.write(chunk)
    finally:
        view_process.stdin.close()
        filter_process.stdin.close()
        bwa_process.stdout.close()
    statuses = {
        "bwa mem": bwa_process.wait(),
        "samtools view": view_process.wait(),
        "samtools fastq": fastq_process.wait(),
        "read recruitment": filter_process.wait(),
    }
    failures = {name: status for name, status in statuses.items() if status}
    if failures:
        raise ReassemblyError(f"Read-recruitment pipeline failed: {failures}")
    (reads_output / ".recruitment_complete").write_text("complete\n")


def run_spades(
    read_file: Path,
    output: Path,
    args: argparse.Namespace,
    cpus: int,
    memory_gb: int,
) -> None:
    bin_name = read_file.name.removesuffix("_1.fastq")
    assembly_dir = output / "reassemblies" / bin_name
    scaffolds = assembly_dir / "scaffolds.fasta"
    if scaffolds.is_file() and scaffolds.stat().st_size:
        print(f"[INFO] Reusing completed SPAdes output: {scaffolds}")
        return
    temp_root = Path(args.tmp_dir).expanduser().resolve()
    temporary = Path(tempfile.mkdtemp(prefix=f"metahict_spades_{bin_name}_", dir=temp_root))
    command: list[object] = [
        required_program("spades.py"),
        "-t",
        cpus,
        "-m",
        memory_gb,
        "--tmp",
        temporary,
    ]
    if args.spades_mode == "careful":
        command.append("--careful")
    if args.spades_phred_offset:
        command.extend(["--phred-offset", args.spades_phred_offset])
    command.extend(shlex.split(args.spades_extra_args))
    command.extend(
        [
            "--untrusted-contigs", output / "original_bins" / f"{bin_name.rsplit('.', 1)[0]}.fa",
            "-1", read_file,
            "-2", read_file.with_name(read_file.name.replace("_1.fastq", "_2.fastq")),
            "-o", assembly_dir,
        ]
    )
    try:
        run(command)
    except subprocess.CalledProcessError:
        print(f"[ERROR] SPAdes temporary directory retained: {temporary}", file=sys.stderr)
        raise
    else:
        if args.keep_temp:
            print(f"[INFO] SPAdes temporary directory retained: {temporary}")
        else:
            shutil.rmtree(temporary, ignore_errors=True)


def parallel_spades_resources(
    jobs: int, threads: int, memory_gb: int
) -> tuple[int, int, int]:
    if jobs < 1 or threads < 1 or memory_gb < 1:
        raise ReassemblyError("SPAdes jobs, threads, and memory must be positive")
    workers = min(jobs, threads, max(1, memory_gb // 8))
    cpus_per_job = max(1, threads // workers)
    memory_per_job = max(1, memory_gb // workers)
    return workers, cpus_per_job, memory_per_job


def run_all_spades(output: Path, args: argparse.Namespace) -> None:
    read_files = sorted((output / "reads_for_reassembly").glob("*_1.fastq"))
    if not read_files:
        raise ReassemblyError("Read recruitment produced no per-bin FASTQ files")
    (output / "reassemblies").mkdir(exist_ok=True)
    if args.parallel:
        workers, cpus_per_job, memory_per_job = parallel_spades_resources(
            len(read_files), args.threads, args.memory
        )
        print(
            "[INFO] Parallel SPAdes allocation: "
            f"{workers} jobs, {cpus_per_job} CPUs/job, {memory_per_job} GB/job"
        )
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    run_spades,
                    path,
                    output,
                    args,
                    cpus_per_job,
                    memory_per_job,
                )
                for path in read_files
            ]
            for future in as_completed(futures):
                future.result()
    else:
        for path in read_files:
            run_spades(path, output, args, args.threads, args.memory)


def finalize_reassemblies(output: Path, minimum_length: int) -> Path:
    destination = output / "reassembled_bins"
    replace_path(destination)
    destination.mkdir()
    for assembly_dir in sorted((output / "reassemblies").iterdir()):
        scaffolds = assembly_dir / "scaffolds.fasta"
        if scaffolds.is_file() and scaffolds.stat().st_size:
            write_minimum_length_fasta(scaffolds, destination / f"{assembly_dir.name}.fa", minimum_length)
    if not any(destination.iterdir()):
        raise ReassemblyError("None of the bins were successfully reassembled")
    return destination


def assemble_residual_reads(output: Path, args: argparse.Namespace) -> Path:
    residual = output / "residual_contigs.fa"
    replace_path(residual)
    forward = output / "unmapped_shotgun_1.fastq"
    reverse = output / "unmapped_shotgun_2.fastq"
    if args.skip_residual_assembly or not forward.is_file() or not reverse.is_file() or not forward.stat().st_size or not reverse.stat().st_size:
        residual.touch()
        return residual
    assembly_dir = output / "residual_assembly"
    replace_path(assembly_dir)
    temporary = Path(tempfile.mkdtemp(prefix="metahict_megahit_residual_", dir=args.tmp_dir))
    try:
        run(
            [
                required_program("megahit"), "-1", forward, "-2", reverse,
                "-o", assembly_dir, "--tmp-dir", temporary,
                "--num-cpu-threads", args.threads, "--min-contig-len", args.minimum_length,
                "--k-min", 21, "--k-max", 141, "--k-step", 12, "--merge-level", "20,0.95",
            ]
        )
    finally:
        if args.keep_temp:
            print(f"[INFO] MEGAHIT temporary directory retained: {temporary}")
        else:
            shutil.rmtree(temporary, ignore_errors=True)
    final_contigs = assembly_dir / "final.contigs.fa"
    if final_contigs.is_file() and final_contigs.stat().st_size:
        shutil.copy2(final_contigs, residual)
    else:
        residual.touch()
    return residual


def checkm2_executable() -> str:
    executable = os.environ.get("METAHICT_CHECKM2_BIN") or shutil.which("checkm2")
    if not executable:
        raise ReassemblyError("CheckM2 was not found; run ./metahict install")
    return executable


def checkm2_environment(executable: str) -> dict[str, str]:
    """Expose CheckM2 companion programs from the executable's environment."""
    environment = os.environ.copy()
    executable_dir = str(Path(executable).expanduser().absolute().parent)
    path_entries = [
        entry
        for entry in environment.get("PATH", "").split(os.pathsep)
        if entry and entry != executable_dir
    ]
    environment["PATH"] = os.pathsep.join([executable_dir, *path_entries])
    if not shutil.which("diamond", path=environment["PATH"]):
        raise ReassemblyError(
            f"DIAMOND was not found in the CheckM2 environment: {executable_dir}"
        )
    return environment


def run_checkm2(input_dir: Path, output_dir: Path, args: argparse.Namespace) -> Path:
    replace_path(output_dir)
    temporary = Path(tempfile.mkdtemp(prefix="metahict_checkm2_", dir=args.tmp_dir))
    executable = checkm2_executable()
    environment = checkm2_environment(executable)
    environment["TMPDIR"] = str(temporary)
    if args.checkm2_db:
        environment["CHECKM2DB"] = str(Path(args.checkm2_db).expanduser().resolve())
    try:
        run(
            [
                executable, "predict", "-i", input_dir, "-o", output_dir,
                "-x", "fa", "-t", args.threads, "--tmpdir", temporary,
            ],
            env=environment,
        )
    finally:
        if args.keep_temp:
            print(f"[INFO] CheckM2 temporary directory retained: {temporary}")
        else:
            shutil.rmtree(temporary, ignore_errors=True)
    report = output_dir / "quality_report.tsv"
    if not report.is_file() or not report.stat().st_size:
        raise ReassemblyError(f"CheckM2 did not create a quality report: {report}")
    return report


def summarize_checkm2(report: Path, output: Path) -> None:
    result = run([sys.executable, HELPERS / "summarize_checkm2.py", report], capture=True)
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        raise ReassemblyError(f"CheckM2 summarizer produced no output for {report}")
    output.write_text("\n".join([lines[0], *sorted(lines[1:])]) + "\n")


def select_best_bins(output: Path, args: argparse.Namespace) -> None:
    candidates = output / "reassembled_bins"
    for original in sorted((output / "original_bins").glob("*.fa")):
        shutil.copy2(original, candidates / f"{original.stem}.orig.fa")
    first_report = run_checkm2(candidates, output / "reassembled_bins.checkm2", args)
    first_stats = output / "candidate_bins_quality.tsv"
    summarize_checkm2(first_report, first_stats)
    result = run(
        [
            sys.executable, HELPERS / "choose_best_bin.py", first_stats,
            args.completeness, args.contamination, args.contamination_penalty,
        ],
        capture=True,
    )
    selected = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    best = output / "reassembled_best_bins"
    replace_path(best)
    best.mkdir()
    name_map = output / "reassembled_bin_name_map.tsv"
    counts = {"orig": 0, "strict": 0, "permissive": 0}
    with name_map.open("w") as mapping:
        mapping.write("final_bin\tselected_candidate\tselection_type\n")
        for candidate in selected:
            source = candidates / f"{candidate}.fa"
            if not source.is_file():
                raise ReassemblyError(f"Selected best-bin candidate is missing: {source}")
            selection_type = candidate.rsplit(".", 1)[-1]
            base_bin = candidate.rsplit(".", 1)[0]
            final_bin = base_bin.replace(".", "")
            shutil.copy2(source, best / f"{final_bin}.fa")
            mapping.write(f"{final_bin}\t{candidate}\t{selection_type}\n")
            counts[selection_type] = counts.get(selection_type, 0) + 1
    if not any(best.iterdir()):
        raise ReassemblyError("No bins passed the final quality thresholds")
    print(
        "[INFO] Best-bin selection: "
        f"{counts['strict']} strict, {counts['permissive']} permissive, {counts['orig']} original"
    )

    work = output / "work_files"
    replace_path(work)
    work.mkdir()
    candidates.replace(work / "reassembled_bins")
    first_stats.replace(work / "candidate_bins_quality.tsv")
    shutil.copy2(name_map, work / name_map.name)
    for name in ("reads_for_reassembly", "binned_assembly", "reassemblies"):
        path = output / name
        if path.exists():
            path.replace(work / name)
    best.replace(output / "reassembled_bins")

    final_report = run_checkm2(output / "reassembled_bins", output / "reassembled_bins.checkm2", args)
    final_stats = output / "reassembled_bins_quality.tsv"
    summarize_checkm2(final_report, final_stats)
    run([sys.executable, HELPERS / "plot_checkm2_results.py", final_report, output / "reassembled_bins"])
    old_header, *old_records = (work / "candidate_bins_quality.tsv").read_text().splitlines()
    quality = output / "quality"
    quality.mkdir(exist_ok=True)
    original_stats = quality / "original_bins_quality.tsv"
    original_stats.write_text(old_header + "\n" + "\n".join(line for line in old_records if "orig" in line) + "\n")
    shutil.copy2(final_report, quality / "checkm2_quality_report.tsv")
    run(
        [
            sys.executable, HELPERS / "plot_reassembly.py", output,
            args.completeness, args.contamination, final_stats, original_stats,
        ]
    )


def write_combined_contigs(bins_directory: Path, residual: Path, output: Path) -> None:
    seen: dict[str, int] = {}
    with output.open("w") as handle:
        for fasta in sorted(bins_directory.iterdir()):
            if fasta.suffix not in {".fa", ".fasta", ".fna"}:
                continue
            bin_name = fasta.stem
            for old_id, sequence in read_fasta(fasta):
                base_id = f"{bin_name}|{old_id}"
                seen[base_id] = seen.get(base_id, 0) + 1
                identifier = base_id if seen[base_id] == 1 else f"{base_id}|copy{seen[base_id]}"
                handle.write(f">{identifier}\n")
                handle.write("\n".join(sequence) + "\n")
        if residual.is_file():
            for old_id, sequence in read_fasta(residual):
                base_id = f"residual|{old_id}"
                seen[base_id] = seen.get(base_id, 0) + 1
                identifier = base_id if seen[base_id] == 1 else f"{base_id}|copy{seen[base_id]}"
                handle.write(f">{identifier}\n")
                handle.write("\n".join(sequence) + "\n")
    if not output.stat().st_size:
        raise ReassemblyError("Combined contigs FASTA is empty")


DEFAULT_THREADS = 16
DEFAULT_MEMORY_GB = 51


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Reassemble METAHICT bins")
    parser.add_argument("-b", "--bins", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-1", dest="forward_reads", required=True)
    parser.add_argument("-2", dest="reverse_reads", required=True)
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        default=DEFAULT_MEMORY_GB,
        help="Total SPAdes memory in GB (default: 51; 80%% of workflow allocation)",
    )
    parser.add_argument("-c", "--completeness", type=float, default=50)
    parser.add_argument("-x", "--contamination", type=float, default=10)
    parser.add_argument("-l", "--minimum-length", type=int, default=500)
    parser.add_argument("--strict-cut-off", type=int, default=2)
    parser.add_argument("--permissive-cut-off", type=int, default=5)
    parser.add_argument("--contamination-penalty", type=float, default=5)
    parser.add_argument("--skip-checkm", dest="run_checkm", action="store_false")
    parser.add_argument("--parallel", action="store_true")
    parser.add_argument("--checkm2_db")
    parser.add_argument("--tmp-dir", default=os.environ.get("METAHICT_TMP_ROOT", os.environ.get("TMPDIR", "/tmp")))
    parser.add_argument("--spades-mode", choices=("careful", "none"), default="careful")
    parser.add_argument("--spades-phred-offset")
    parser.add_argument("--spades-extra-args", default="")
    parser.add_argument("--skip-residual-assembly", action="store_true")
    parser.add_argument("--keep-temp", action="store_true")
    parser.set_defaults(run_checkm=True)
    return parser


def command_main(args: argparse.Namespace) -> Path:
    bins = Path(args.bins).expanduser().resolve()
    forward_reads = Path(args.forward_reads).expanduser().resolve()
    reverse_reads = Path(args.reverse_reads).expanduser().resolve()
    if not bins.is_dir():
        raise ReassemblyError(f"Bin directory not found: {bins}")
    for path in (forward_reads, reverse_reads):
        if not path.is_file():
            raise ReassemblyError(f"Read file not found: {path}")
    temp_root = Path(args.tmp_dir).expanduser().resolve()
    temp_root.mkdir(parents=True, exist_ok=True)
    args.tmp_dir = str(temp_root)
    output = Path(args.output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    original_bins = output / "original_bins"
    replace_path(original_bins)
    shutil.copytree(bins, original_bins)
    binned_assembly = output / "binned_assembly"
    binned_assembly.mkdir(exist_ok=True)
    assembly = binned_assembly / "assembly.fa"
    if concatenate_fastas(original_bins, assembly) == 0 or not assembly.stat().st_size:
        raise ReassemblyError(f"No input bin FASTAs were found in {original_bins}")
    recruitment_marker = output / "reads_for_reassembly" / ".recruitment_complete"
    if not recruitment_marker.is_file():
        recruit_reads(
            assembly, forward_reads, reverse_reads, original_bins, output,
            args.threads, args.strict_cut_off, args.permissive_cut_off,
        )
    else:
        print(f"[INFO] Reusing completed recruited reads: {recruitment_marker.parent}")
    run_all_spades(output, args)
    reassembled_bins = finalize_reassemblies(output, args.minimum_length)
    residual = assemble_residual_reads(output, args)
    if args.run_checkm:
        select_best_bins(output, args)
        reassembled_bins = output / "reassembled_bins"
    write_combined_contigs(reassembled_bins, residual, output / "combined_contigs.fa")
    print(f"[PASS] Reassembled bins: {reassembled_bins}")
    return reassembled_bins


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        command_main(build_parser().parse_args(argv))
    except (ReassemblyError, OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
