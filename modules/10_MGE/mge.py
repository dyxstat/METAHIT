#!/usr/bin/env python3
"""Detect mobile genetic elements, find circular contigs, and report MGE–host pairs."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import Iterable, Optional, Sequence

np = None
pd = None
SeqIO = None
scipy_load_npz = None


class MgeError(RuntimeError):
    """Expected user-facing MGE failure."""


def load_scientific_dependencies() -> None:
    """Import analysis dependencies only after argument parsing.

    Keeping these imports lazy lets ``mge.py --help`` work before the locked
    runtime is installed while still failing clearly for a real run.
    """
    global np, pd, SeqIO, scipy_load_npz
    try:
        import numpy as numpy_module
        import pandas as pandas_module
        from Bio import SeqIO as seqio_module
    except ModuleNotFoundError as error:
        raise MgeError(
            f"Missing scientific Python dependency '{error.name}'. Run ./metahict install."
        ) from error
    try:
        from scipy.sparse import load_npz as sparse_loader
    except ModuleNotFoundError:
        sparse_loader = None
    np = numpy_module
    pd = pandas_module
    SeqIO = seqio_module
    scipy_load_npz = sparse_loader


class SimpleCOO:
    def __init__(self, row, col, data, shape):
        self.row = row
        self.col = col
        self.data = data
        self.shape = tuple(int(value) for value in shape)


def run(command: Sequence[object], *, env: Optional[dict[str, str]] = None) -> None:
    argv = [str(item) for item in command]
    print(f"[RUN] {shlex.join(argv)}", flush=True)
    subprocess.run(argv, check=True, env=env)


def existing_file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise MgeError(f"{label} not found: {path}")
    return path


def executable(name: str, candidates: Iterable[Path]) -> str:
    found = shutil.which(name)
    if found:
        return found
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    raise MgeError(f"Required executable not found: {name}")


def ccfind_runtime_environment(
    project: Path, base_environment: Optional[dict[str, str]] = None
) -> dict[str, str]:
    """Expose ccfind and its locked FASTA-suite companion programs."""
    source_environment = os.environ if base_environment is None else base_environment
    environment = source_environment.copy()
    ccfind_bin = project / "conda_envs" / "ccfind_env" / "bin"
    if ccfind_bin.is_dir():
        environment["PATH"] = os.pathsep.join(
            [str(ccfind_bin), environment.get("PATH", "")]
        )
    if shutil.which("ssearch36", path=environment.get("PATH", "")) is None:
        raise MgeError(
            "ccfind companion executable 'ssearch36' was not found. "
            "Verify conda_envs/ccfind_env with './metahict verify'."
        )
    return environment


def prepare_unique_fasta(source: Path, destination: Path, map_file: Path) -> tuple[int, int]:
    """Copy FASTA while deterministically renaming duplicate primary IDs."""
    seen: dict[str, int] = {}
    used: set[str] = set()
    total = 0
    renamed = 0
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open() as src, destination.open("w") as out, map_file.open("w") as mapping:
        mapping.write("original_id\tunique_id\toccurrence\n")
        for line in src:
            if not line.startswith(">"):
                out.write(line)
                continue
            total += 1
            header = line[1:].rstrip("\n")
            parts = header.split(maxsplit=1)
            original_id = parts[0]
            description = f" {parts[1]}" if len(parts) > 1 else ""
            occurrence = seen.get(original_id, 0) + 1
            seen[original_id] = occurrence
            unique_id = original_id
            if occurrence > 1 or unique_id in used:
                renamed += 1
                suffix = occurrence
                unique_id = f"{original_id}__dup{suffix}"
                while unique_id in used:
                    suffix += 1
                    unique_id = f"{original_id}__dup{suffix}"
            used.add(unique_id)
            mapping.write(f"{original_id}\t{unique_id}\t{occurrence}\n")
            out.write(f">{unique_id}{description}\n")
    return total, renamed


def filter_proviruses(virus_summary: Path, viral_fasta: Path, output: Path) -> None:
    summary = pd.read_csv(virus_summary, sep="\t")
    if "topology" not in summary.columns or "seq_name" not in summary.columns:
        raise MgeError(f"geNomad virus summary lacks topology/seq_name columns: {virus_summary}")
    keep_ids = set(summary.loc[summary["topology"] != "Provirus", "seq_name"].astype(str))
    with output.open("w") as handle:
        for record in SeqIO.parse(str(viral_fasta), "fasta"):
            if record.id in keep_ids:
                SeqIO.write(record, handle, "fasta")


def prepare_ccfind_fasta(source: Path, destination: Path, map_file: Path) -> int:
    destination.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with destination.open("w") as out, map_file.open("w") as mapping:
        mapping.write("ccfind_id\tcontig_id\n")
        for count, record in enumerate(SeqIO.parse(str(source), "fasta"), 1):
            safe_id = f"ccfind_{count:09d}"
            mapping.write(f"{safe_id}\t{record.id}\n")
            record.id = safe_id
            record.name = safe_id
            record.description = safe_id
            SeqIO.write(record, out, "fasta")
    return count


def remap_ccfind_ids(map_file: Path, raw_file: Path, mapped_file: Path) -> None:
    safe_to_contig: dict[str, str] = {}
    with map_file.open() as handle:
        next(handle, None)
        for line in handle:
            safe_id, contig_id = line.rstrip("\n").split("\t", 1)
            safe_to_contig[safe_id] = contig_id
    with raw_file.open() as source, mapped_file.open("w") as output:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            fields[0] = safe_to_contig.get(fields[0], fields[0])
            output.write("\t".join(fields) + "\n")


def load_npz_coo(path: Path):
    if scipy_load_npz is not None:
        return scipy_load_npz(path).tocoo()
    archive = np.load(path, allow_pickle=False)
    fmt = archive["format"]
    if hasattr(fmt, "item"):
        fmt = fmt.item()
    if isinstance(fmt, bytes):
        fmt = fmt.decode()
    shape = tuple(int(value) for value in archive["shape"])
    data = archive["data"]
    if str(fmt) == "coo":
        return SimpleCOO(archive["row"], archive["col"], data, shape)
    if str(fmt) == "csr":
        row = np.repeat(np.arange(shape[0]), np.diff(archive["indptr"]))
        return SimpleCOO(row, archive["indices"], data, shape)
    if str(fmt) == "csc":
        col = np.repeat(np.arange(shape[1]), np.diff(archive["indptr"]))
        return SimpleCOO(archive["indices"], col, data, shape)
    raise MgeError(f"Unsupported sparse matrix format in {path}: {fmt}")


def count_host_unmapped(
    fasta: Path, contig_to_host: dict[str, str]
) -> tuple[int, int, int]:
    identifiers = [record.id for record in SeqIO.parse(str(fasta), "fasta")]
    hosts = sum(identifier in contig_to_host for identifier in identifiers)
    return len(identifiers), hosts, len(identifiers) - hosts


def fasta_lengths(fasta: Path) -> dict[str, int]:
    if not fasta.exists():
        return {}
    return {record.id: len(record.seq) for record in SeqIO.parse(str(fasta), "fasta")}


def fasta_sequences(fasta: Path) -> dict[str, str]:
    if not fasta.exists():
        return {}
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(str(fasta), "fasta")}


def read_ccfind_detected(path: Path) -> set[str]:
    if not path.is_file() or path.stat().st_size == 0:
        return set()
    with path.open() as handle:
        return {
            line.strip().split("\t", 1)[0]
            for line in handle
            if line.strip() and not line.startswith("#")
        }


def count_nonempty_lines(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    with path.open() as handle:
        return sum(bool(line.strip()) and not line.startswith("#") for line in handle)


def read_summary(path: Path) -> pd.DataFrame:
    if path.is_file() and path.stat().st_size:
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def summary_id_column(frame: pd.DataFrame) -> Optional[str]:
    for column in ("seq_name", "contig_id", "id"):
        if column in frame.columns:
            return column
    return str(frame.columns[0]) if len(frame.columns) else None


def topology_map(frame: pd.DataFrame) -> dict[str, str]:
    id_column = summary_id_column(frame)
    if id_column is None or "topology" not in frame.columns:
        return {}
    return dict(zip(frame[id_column].astype(str), frame["topology"].fillna("").astype(str)))


def genomad_circular_evidence(topology: object) -> bool:
    return any(term in str(topology).lower() for term in ("circular", "dtr", "itr"))


def fasta_primary_ids(path: Path) -> list[str]:
    """Read primary FASTA identifiers without requiring scientific packages."""
    identifiers = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            identifier = line[1:].strip().split(maxsplit=1)[0]
            if not identifier:
                raise MgeError(f"Empty FASTA identifier in {path}")
            identifiers.append(identifier)
    return identifiers


def build_contig_to_host(
    fasta: Path, host_dir: Path, output: Path
) -> dict[str, str]:
    """Build and validate explicit FASTA-contig to host-MAG membership."""
    if not host_dir.is_dir():
        raise MgeError(f"Host-genome directory not found: {host_dir}")
    assembly_ids = fasta_primary_ids(fasta)
    if not assembly_ids:
        raise MgeError(f"Metagenome FASTA contains no records: {fasta}")
    identifier_counts: dict[str, int] = {}
    for identifier in assembly_ids:
        identifier_counts[identifier] = identifier_counts.get(identifier, 0) + 1
    duplicate_assembly_ids = sorted(
        identifier for identifier, count in identifier_counts.items() if count > 1
    )
    if duplicate_assembly_ids:
        preview = ", ".join(duplicate_assembly_ids[:5])
        raise MgeError(
            "MGE–host pairing requires unique FASTA identifiers; "
            f"duplicates include: {preview}"
        )
    assembly_id_set = set(assembly_ids)
    fasta_files = sorted(
        path
        for path in host_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna"}
    )
    if not fasta_files:
        raise MgeError(
            f"Host-genome directory contains no .fa, .fasta, or .fna files: {host_dir}"
        )

    contig_to_host: dict[str, str] = {}
    unmatched: list[str] = []
    for fasta in fasta_files:
        host_id = fasta.stem
        host_contigs = fasta_primary_ids(fasta)
        if not host_contigs:
            raise MgeError(f"Host-genome FASTA contains no records: {fasta}")
        for contig_id in host_contigs:
            candidates = [
                candidate
                for candidate in (contig_id, f"{host_id}|{contig_id}")
                if candidate in assembly_id_set
            ]
            if len(candidates) != 1:
                unmatched.append(f"{host_id}:{contig_id}")
                continue
            assembly_id = candidates[0]
            previous = contig_to_host.get(assembly_id)
            if previous is not None and previous != host_id:
                raise MgeError(
                    f"Metagenome contig belongs to multiple hosts: {assembly_id} "
                    f"({previous}, {host_id})"
                )
            contig_to_host[assembly_id] = host_id
    if unmatched:
        preview = ", ".join(unmatched[:10])
        raise MgeError(
            "Host contigs could not be matched uniquely to the metagenome FASTA; "
            f"examples: {preview}"
        )
    if not contig_to_host:
        raise MgeError("No assembly contigs were assigned to host MAGs")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("contig_id\thost_id\n")
        for contig_id in assembly_ids:
            if contig_id in contig_to_host:
                handle.write(f"{contig_id}\t{contig_to_host[contig_id]}\n")
    return contig_to_host


def extract_taxonomy_map(frame: pd.DataFrame) -> dict[str, str]:
    if len(frame.columns) == 0:
        return {}
    id_column = "seq_name" if "seq_name" in frame.columns else str(frame.columns[0])
    taxonomy_columns = [
        column
        for column in frame.columns
        if any(token in column.lower() for token in ("tax", "lineage", "virus_name", "classification"))
    ]
    if not taxonomy_columns:
        return {}
    tax_frame = frame[[id_column] + taxonomy_columns].copy()
    tax_frame["mge_taxonomy"] = tax_frame[taxonomy_columns].fillna("").astype(str).agg(";".join, axis=1)
    return dict(zip(tax_frame[id_column].astype(str), tax_frame["mge_taxonomy"]))


def contact_contigs_for_matrix(matrix_file: Path, combined_fasta: Path) -> tuple[list[str], Path]:
    contig_info = matrix_file.resolve().parent / "contig_info.csv"
    if contig_info.exists():
        frame = pd.read_csv(contig_info)
        if "name" not in frame.columns:
            raise MgeError(f"Contact contig info file lacks a 'name' column: {contig_info}")
        return frame["name"].astype(str).tolist(), contig_info
    identifiers = [record.id for record in SeqIO.parse(str(combined_fasta), "fasta")]
    return identifiers, combined_fasta


def write_reports(
    *,
    outdir: Path,
    virus_summary: Path,
    plasmid_summary: Path,
    filtered_viral: Path,
    plasmid_fasta: Path,
    contact_matrix: Path,
    raw_contact_matrix: Path,
    assembly: Path,
    contig_to_host: dict[str, str],
    pair_filter: str,
    zscore_threshold: float,
    fixed_contact_threshold: float,
    top_percent: float,
    min_raw_contacts: float,
    min_contact_strength: float,
    ccfind_detected_list: Path,
    ccfind_too_short_list: Path,
    ccfind_terminal_fragment_size: int,
    ccfind_min_identity: int,
    ccfind_min_aligned_length: int,
) -> Path:
    summary_file = outdir / "mge_reports" / "MGE_summary.txt"
    virus_frame = pd.read_csv(virus_summary, sep="\t")
    plasmid_frame = read_summary(plasmid_summary)
    virus_ids = set(fasta_lengths(filtered_viral))
    plasmid_ids = set(fasta_lengths(plasmid_fasta))
    virus_topology = topology_map(virus_frame)
    plasmid_topology = topology_map(plasmid_frame)
    ccfind_detected = read_ccfind_detected(ccfind_detected_list)
    assembly_sequences = fasta_sequences(assembly)
    mge_ids = virus_ids | plasmid_ids

    topology_rows = []
    for contig_id in sorted(assembly_sequences):
        genomad: object = "NA"
        if contig_id in virus_ids:
            genomad = int(genomad_circular_evidence(virus_topology.get(contig_id, "")))
        elif contig_id in plasmid_ids:
            genomad = int(genomad_circular_evidence(plasmid_topology.get(contig_id, "")))
        topology_rows.append(
            {
                "contig_id": contig_id,
                "mge": int(contig_id in mge_ids),
                "ccfind": int(contig_id in ccfind_detected),
                "genomad": genomad,
            }
        )
    topology_columns = ["contig_id", "mge", "ccfind", "genomad"]
    topology_frame = pd.DataFrame(topology_rows).reindex(columns=topology_columns)
    topology_frame.to_csv(outdir / "sequence_topology.tsv", sep="\t", index=False)

    mge_type_map = {contig: "virus" for contig in virus_ids}
    mge_type_map.update({contig: "plasmid" for contig in plasmid_ids if contig not in mge_type_map})
    mge_contigs = set(mge_type_map)
    taxonomy_map = extract_taxonomy_map(virus_frame)
    taxonomy_map.update(extract_taxonomy_map(plasmid_frame))
    normalized_pairs: dict[tuple[str, str], dict[str, object]] = {}
    raw_support: dict[tuple[str, str], float] = {}
    normalized_matrix = load_npz_coo(contact_matrix)
    raw_matrix = load_npz_coo(raw_contact_matrix)
    contigs, contig_order_source = contact_contigs_for_matrix(contact_matrix, assembly)
    raw_contigs, raw_contig_order_source = contact_contigs_for_matrix(raw_contact_matrix, assembly)
    if raw_contigs != contigs:
        raise MgeError("Raw and normalized contact matrices use different contig orders")
    if len(contigs) != len(set(contigs)):
        raise MgeError("Contact-matrix contig order contains duplicate identifiers")
    unknown_contact_contigs = sorted(set(contigs) - set(assembly_sequences))
    if unknown_contact_contigs:
        preview = ", ".join(unknown_contact_contigs[:10])
        raise MgeError(
            "Contact matrix contains identifiers that are absent from the "
            f"selected FASTA; examples: {preview}"
        )
    if normalized_matrix.shape != (len(contigs), len(contigs)):
        raise MgeError(
            f"Contact matrix shape {normalized_matrix.shape} does not match contig count {len(contigs)}"
        )
    if raw_matrix.shape != (len(raw_contigs), len(raw_contigs)):
        raise MgeError(
            f"Raw matrix shape {raw_matrix.shape} does not match contig count {len(raw_contigs)}"
        )
    host_contigs = {
        contig for contig in contigs
        if contig in contig_to_host and contig not in mge_contigs
    }

    def candidate_edges(matrix, threshold: float):
        for row, column, value in zip(matrix.row, matrix.col, matrix.data):
            if row >= column or float(value) <= threshold:
                continue
            first, second = contigs[row], contigs[column]
            if first in mge_contigs and second in host_contigs:
                yield first, second, float(value)
            elif second in mge_contigs and first in host_contigs:
                yield second, first, float(value)

    for mge_contig, host_contig, value in candidate_edges(
        normalized_matrix, min_contact_strength
    ):
        host_mag = contig_to_host[host_contig]
        key = (mge_contig, host_mag)
        pair = normalized_pairs.setdefault(
            key,
            {
                "mge_contig": mge_contig,
                "mge_type": mge_type_map.get(mge_contig, ""),
                "host_mag": host_mag,
                "_host_contigs": set(),
                "normalized_contact_score": 0.0,
                "max_contig_normalized_contact_score": 0.0,
                "mge_taxonomy": taxonomy_map.get(mge_contig, ""),
            },
        )
        pair["_host_contigs"].add(host_contig)  # type: ignore[union-attr]
        pair["normalized_contact_score"] = (
            float(pair["normalized_contact_score"]) + value
        )
        pair["max_contig_normalized_contact_score"] = max(
            float(pair["max_contig_normalized_contact_score"]), value
        )

    for mge_contig, host_contig, value in candidate_edges(raw_matrix, 0):
        key = (mge_contig, contig_to_host[host_contig])
        raw_support[key] = raw_support.get(key, 0.0) + value
    for key, pair in normalized_pairs.items():
        pair["host_contig_count"] = len(pair.pop("_host_contigs"))  # type: ignore[arg-type]
        pair["raw_contact_support"] = raw_support.get(key, 0.0)

    columns = [
        "mge_contig", "mge_type", "host_mag", "raw_contact_support",
        "normalized_contact_score", "max_contig_normalized_contact_score",
        "mge_taxonomy", "host_contig_count", "z_score", "passes_raw_contact_filter",
        "passes_association_filter", "passes_zscore_threshold",
        "association_score_source", "association_filter",
    ]
    pairs = pd.DataFrame(normalized_pairs.values())
    if len(pairs):
        pairs["z_score"] = 0.0
        for indexes in pairs.groupby("mge_type").groups.values():
            scores = pairs.loc[indexes, "normalized_contact_score"]
            standard_deviation = scores.std(ddof=0)
            if standard_deviation != 0 and not pd.isna(standard_deviation):
                pairs.loc[indexes, "z_score"] = (scores - scores.mean()) / standard_deviation
        pairs["passes_raw_contact_filter"] = pairs["raw_contact_support"] >= min_raw_contacts
        if pair_filter == "zscore":
            pair_pass = pairs["z_score"] > zscore_threshold
        elif pair_filter == "fixed":
            pair_pass = pairs["normalized_contact_score"] >= fixed_contact_threshold
        elif pair_filter == "percentage":
            keep = max(1, math.ceil(len(pairs) * top_percent / 100))
            pair_pass = pairs["normalized_contact_score"].rank(method="first", ascending=False) <= keep
        else:
            pair_pass = True
        # Preserve the established result-table column names so existing
        # scientific comparisons and downstream readers remain compatible.
        pairs["passes_association_filter"] = pairs["passes_raw_contact_filter"] & pair_pass
        pairs["passes_zscore_threshold"] = pairs["z_score"] > zscore_threshold
        pairs["association_score_source"] = "normalized"
        pairs["association_filter"] = pair_filter
        pairs = pairs.sort_values(
            ["passes_association_filter", "z_score", "normalized_contact_score"],
            ascending=[False, False, False],
        ).reindex(columns=columns)
        filtered = pairs[pairs["passes_association_filter"]]
    else:
        filtered = pd.DataFrame(columns=columns)
    filtered_file = outdir / f"candidate_mge_host_pairs_{pair_filter}_filtered.tsv"
    filtered.to_csv(filtered_file, sep="\t", index=False)

    virus_total, virus_hosts, virus_unmapped = count_host_unmapped(
        filtered_viral, contig_to_host
    )
    plasmid_total, plasmid_hosts, plasmid_unmapped = count_host_unmapped(
        plasmid_fasta, contig_to_host
    )
    provirus_count = int((virus_frame["topology"] == "Provirus").sum())
    circular_viruses = int(
        ((topology_frame["contig_id"].isin(virus_ids)) & ((topology_frame["ccfind"] == 1) | (topology_frame["genomad"] == 1))).sum()
    )
    circular_plasmids = int(
        ((topology_frame["contig_id"].isin(plasmid_ids)) & ((topology_frame["ccfind"] == 1) | (topology_frame["genomad"] == 1))).sum()
    )
    with summary_file.open("w") as handle:
        handle.write("=== STEP 1.5 Report (After Provirus Removal) ===\n")
        handle.write(f"Virus contigs (after provirus removal): {virus_total}\n  assigned to hosts: {virus_hosts}\n  unassigned: {virus_unmapped}\n")
        handle.write(f"Provirus removed: {provirus_count}\n")
        handle.write(f"Plasmid contigs (raw): {plasmid_total}\n  assigned to hosts: {plasmid_hosts}\n  unassigned: {plasmid_unmapped}\n\n")
        handle.write("=== STEP 2 Report (Detected MGE Candidates) ===\n")
        handle.write(f"Virus candidates: {len(virus_ids)}\n")
        handle.write(f"Plasmid candidates: {len(plasmid_ids)}\n\n")
        handle.write("=== STEP 2.5 Report (Sequence Topology) ===\n")
        handle.write(
            "ccfind parameters: "
            f"terminal_fragment_size={ccfind_terminal_fragment_size}, "
            f"min_percent_identity={ccfind_min_identity}, min_aligned_length={ccfind_min_aligned_length}\n"
        )
        handle.write(f"ccfind too-short contigs excluded: {count_nonempty_lines(ccfind_too_short_list)}\n")
        handle.write(f"Assembly contigs evaluated: {len(topology_rows)}\n")
        handle.write(f"Circular viral candidates by geNomad or ccfind: {circular_viruses}\n")
        handle.write(f"Circular plasmid candidates by geNomad or ccfind: {circular_plasmids}\n")
        handle.write(f"Circular MGE candidates by geNomad topology: {int((topology_frame['genomad'] == 1).sum())}\n")
        handle.write(f"Circular MGE candidates by ccfind: {int(((topology_frame['mge'] == 1) & (topology_frame['ccfind'] == 1)).sum())}\n")
        handle.write(f"ccfind circular non-MGE contigs: {int(((topology_frame['mge'] == 0) & (topology_frame['ccfind'] == 1)).sum())}\n")
        handle.write("Sequence topology report saved to: sequence_topology.tsv\n\n")
        handle.write("=== STEP 3 Report (Candidate MGE–host Pairs: Contact Inputs) ===\n")
        handle.write(
            "Score source for MGE–host pairs: normalized\n"
            f"Normalized contact matrix for MGE–host pairs: {contact_matrix}\n"
        )
        handle.write(f"Raw contact matrix for MGE–host pairs: {raw_contact_matrix}\n")
        handle.write(
            f"Contact contig order for MGE–host pairs: {contig_order_source}\n"
            f"Raw-contact contig order for MGE–host pairs: {raw_contig_order_source}\n"
        )
        handle.write(
            f"Filter for MGE–host pairs: {pair_filter}\n"
            f"Minimum normalized contact strength: {min_contact_strength}\n"
        )
        handle.write(
            f"Minimum raw Hi-C support: {min_raw_contacts}\n"
            f"Viral contigs: {len(virus_ids)}\n"
            f"Plasmid contigs: {len(plasmid_ids)}\n\n"
        )
        handle.write("=== STEP 4 Report (Candidate MGE–host Pairs Table) ===\n")
        if len(pairs):
            handle.write(f"Score source for MGE–host pairs: normalized\nFilter for MGE–host pairs: {pair_filter}\n")
            handle.write(f"Z-score threshold: {zscore_threshold}\nFixed contact threshold: {fixed_contact_threshold}\n")
            handle.write(f"Top percent: {top_percent}\nMinimum raw Hi-C support: {min_raw_contacts}\n")
            handle.write(f"Scored candidate MGE–host pairs before filtering: {len(pairs)}\n")
            handle.write(f"Filtered candidate MGE–host pairs: {len(filtered)}\n")
            handle.write(f"Unique MGE contigs after filtering: {filtered['mge_contig'].nunique()}\n")
            for mge_type, group in filtered.groupby("mge_type"):
                handle.write(f"Unique {mge_type} MGE contigs after filtering: {group['mge_contig'].nunique()}\n")
            handle.write(f"Unique host MAGs after filtering: {filtered['host_mag'].nunique()}\n")
        else:
            handle.write("No candidate MGE–host pairs found.\n")
        handle.write(f"Filtered table of MGE–host pairs saved to: {filtered_file.name}\n")
    return filtered_file


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module mge --help",
    )
    parser.add_argument("-p", "--metahict-path", "--project-path", dest="project_path", required=True, help="METAHICT repository root")
    parser.add_argument(
        "--fasta", required=True,
        help="Metagenome FASTA containing the MGE and host contigs",
    )
    parser.add_argument(
        "--host-dir", required=True,
        help="Directory containing host-genome FASTA files",
    )
    parser.add_argument("--contact", required=True, help="Normalized contact matrix for MGE–host pairing")
    parser.add_argument("--raw-contact", required=True, help="Raw contact matrix for MGE–host pairing")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Threads passed to geNomad and ccfind")
    parser.add_argument("--genomad-db", "--genomad_db", dest="genomad_db", help="geNomad reference database directory")
    parser.add_argument("--genomad-splits", type=int, default=8, help="MMseqs2 database split count")
    parser.add_argument("--genomad-sensitivity", type=float, default=4.2, help="geNomad MMseqs2 sensitivity")
    parser.add_argument("--genomad-cleanup", action=argparse.BooleanOptionalAction, default=True, help="Remove geNomad intermediate files")
    parser.add_argument("--genomad-restart", action="store_true", help="Overwrite existing geNomad intermediate files")
    parser.add_argument("--genomad-preset", choices=("default", "conservative", "relaxed"), default="default", help="geNomad filtering preset")
    parser.add_argument("--genomad-min-score", type=float, default=0.7, help="Minimum virus/plasmid score with the default preset")
    parser.add_argument("--genomad-max-fdr", type=float, default=0.1, help="Maximum false-discovery rate with the default preset")
    parser.add_argument("--genomad-extra-args", default="", help="Additional native geNomad arguments")
    parser.add_argument("--pair-filter", choices=("zscore", "fixed", "percentage", "raw-support-only"), default="zscore", help="Filtering method for MGE–host pairs")
    parser.add_argument("--zscore-threshold", type=float, default=0.5, help="Minimum Z-score for MGE–host pairs")
    parser.add_argument("--fixed-contact-threshold", type=float, default=0, help="Minimum normalized contact in fixed mode")
    parser.add_argument("--top-percent", type=float, default=50, help="Strongest contact percentage retained")
    parser.add_argument("--min-raw-contacts", type=float, default=2, help="Minimum raw Hi-C pair support")
    parser.add_argument("--ccfind-terminal-fragment-size", type=int, default=500, help="ccfind terminal-fragment length")
    parser.add_argument("--ccfind-min-identity", type=int, default=94, help="ccfind minimum terminal-alignment identity")
    parser.add_argument("--ccfind-min-aligned-length", type=int, default=50, help="ccfind minimum aligned length")
    parser.add_argument("--min-contact-strength", type=float, default=0, help="Final minimum normalized contact strength")
    parser.add_argument("--tmp-dir", default=os.environ.get("METAHICT_TMP_ROOT", os.environ.get("TMPDIR", "/tmp")), help="Temporary directory root")
    return parser


def command_main(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise MgeError("--threads must be at least 1")
    if args.genomad_splits < 1 or args.genomad_sensitivity <= 0:
        raise MgeError("geNomad splits must be positive and sensitivity must be greater than 0")
    if not 0 <= args.genomad_min_score <= 1 or not 0 <= args.genomad_max_fdr <= 1:
        raise MgeError("geNomad minimum score and maximum FDR must be between 0 and 1")
    if not 0 < args.top_percent <= 100:
        raise MgeError("--top-percent must be > 0 and <= 100")
    if min(
        args.min_raw_contacts,
        args.min_contact_strength,
        args.ccfind_terminal_fragment_size,
        args.ccfind_min_identity,
        args.ccfind_min_aligned_length,
    ) < 0:
        raise MgeError("Thresholds for MGE–host pairs and ccfind must be non-negative")
    if args.ccfind_min_identity > 100:
        raise MgeError("--ccfind-min-identity must be at most 100")
    load_scientific_dependencies()
    project = Path(args.project_path).expanduser().resolve()
    assembly_source = existing_file(args.fasta, "Metagenome FASTA")
    contact_matrix = existing_file(args.contact, "Normalized contact matrix")
    raw_contact_matrix = existing_file(args.raw_contact, "Raw contact matrix")
    outdir = Path(args.outdir).expanduser().resolve()
    temp_root = Path(args.tmp_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    temp_root.mkdir(parents=True, exist_ok=True)
    (outdir / "genomad_output").mkdir(exist_ok=True)
    (outdir / "mge_reports").mkdir(exist_ok=True)
    genomad_database = Path(args.genomad_db).expanduser().resolve() if args.genomad_db else project / "databases" / "genomad_db" / "genomad_db"
    if not genomad_database.is_dir():
        raise MgeError(f"geNomad database missing: {genomad_database}")
    genomad = executable("genomad", [project / "conda_envs" / "genomad" / "bin" / "genomad"])
    ccfind = executable(
        "ccfind",
        [project / "external" / "bin" / "ccfind", project / "conda_envs" / "ccfind_env" / "bin" / "ccfind"],
    )
    assembly_dir = outdir / "input_assembly"
    assembly = assembly_dir / "assembly.unique.fa"
    id_map = assembly_dir / "assembly.id_map.tsv"
    total, renamed = prepare_unique_fasta(assembly_source, assembly, id_map)
    print(f"[INFO] Assembly FASTA records: {total}; duplicate IDs renamed: {renamed}")
    if renamed:
        raise MgeError(
            "MGE–host pairing requires unique FASTA identifiers so "
            "contact and host mappings remain unambiguous"
        )
    contig_to_host = build_contig_to_host(
        assembly,
        Path(args.host_dir).expanduser().resolve(),
        assembly_dir / "contig_to_host.tsv",
    )

    genomad_command: list[object] = [
        genomad, "end-to-end", "--splits", args.genomad_splits,
        "--sensitivity", args.genomad_sensitivity, "--threads", args.threads,
    ]
    if args.genomad_cleanup:
        genomad_command.append("--cleanup")
    if args.genomad_restart:
        genomad_command.append("--restart")
    if args.genomad_preset == "conservative":
        genomad_command.append("--conservative")
    elif args.genomad_preset == "relaxed":
        genomad_command.append("--relaxed")
    else:
        genomad_command.extend(["--min-score", args.genomad_min_score, "--max-fdr", args.genomad_max_fdr])
    genomad_command.extend(shlex.split(args.genomad_extra_args))
    genomad_command.extend([assembly, outdir / "genomad_output", genomad_database])
    environment = os.environ.copy()
    environment["TMPDIR"] = str(temp_root)
    run(genomad_command, env=environment)

    assembly_base = assembly.stem
    summary_dir = outdir / "genomad_output" / f"{assembly_base}_summary"
    virus_summary = existing_file(str(summary_dir / f"{assembly_base}_virus_summary.tsv"), "geNomad virus summary")
    plasmid_summary = existing_file(str(summary_dir / f"{assembly_base}_plasmid_summary.tsv"), "geNomad plasmid summary")
    viral_fasta = existing_file(str(summary_dir / f"{assembly_base}_virus.fna"), "geNomad virus FASTA")
    plasmid_fasta = existing_file(str(summary_dir / f"{assembly_base}_plasmid.fna"), "geNomad plasmid FASTA")
    filtered_viral = outdir / "mge_reports" / "virus_no_provirus.fna"
    filter_proviruses(virus_summary, viral_fasta, filtered_viral)

    ccfind_output = outdir / "ccfind_output"
    if ccfind_output.exists():
        shutil.rmtree(ccfind_output)
    ccfind_input_dir = outdir / "ccfind_input"
    ccfind_input = ccfind_input_dir / "assembly.ccfind_safe.fa"
    ccfind_map = ccfind_input_dir / "ccfind_id_map.tsv"
    prepare_ccfind_fasta(assembly, ccfind_input, ccfind_map)
    ccfind_command: list[object] = [
        ccfind, ccfind_input, ccfind_output,
        "--terminal-fragment-size", args.ccfind_terminal_fragment_size,
        "--min-percent-identity", args.ccfind_min_identity,
        "--min-aligned-length", args.ccfind_min_aligned_length,
    ]
    if args.threads > 1:
        ccfind_command.extend(["--ncpus", args.threads])
    run(ccfind_command, env=ccfind_runtime_environment(project, environment))
    result_dir = ccfind_output / "result"
    result_dir.mkdir(parents=True, exist_ok=True)
    raw_detected = result_dir / "circ.detected.list"
    too_short = result_dir / "too_short_seq.list"
    raw_detected.touch(exist_ok=True)
    too_short.touch(exist_ok=True)
    mapped_detected = result_dir / "circ.detected.mapped.list"
    remap_ccfind_ids(ccfind_map, raw_detected, mapped_detected)

    filtered_file = write_reports(
        outdir=outdir,
        virus_summary=virus_summary,
        plasmid_summary=plasmid_summary,
        filtered_viral=filtered_viral,
        plasmid_fasta=plasmid_fasta,
        contact_matrix=contact_matrix,
        raw_contact_matrix=raw_contact_matrix,
        assembly=assembly,
        contig_to_host=contig_to_host,
        pair_filter=args.pair_filter,
        zscore_threshold=args.zscore_threshold,
        fixed_contact_threshold=args.fixed_contact_threshold,
        top_percent=args.top_percent,
        min_raw_contacts=args.min_raw_contacts,
        min_contact_strength=args.min_contact_strength,
        ccfind_detected_list=mapped_detected,
        ccfind_too_short_list=too_short,
        ccfind_terminal_fragment_size=args.ccfind_terminal_fragment_size,
        ccfind_min_identity=args.ccfind_min_identity,
        ccfind_min_aligned_length=args.ccfind_min_aligned_length,
    )
    # The remapped ccfind result is retained under ccfind_output/. Its rewritten
    # FASTA and identifier map are reproducible process inputs, not final results.
    shutil.rmtree(ccfind_input_dir)
    print(f"[PASS] MGE analysis complete: {filtered_file}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        command_main(build_parser().parse_args(argv))
    except (MgeError, OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
