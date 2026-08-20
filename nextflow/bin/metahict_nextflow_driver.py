#!/usr/bin/env python3
import argparse
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


MODULE_ORDER = [
    "preprocessing_sg",
    "preprocessing_hic",
    "assembly",
    "alignment",
    "coverage",
    "contact",
    "binning",
    "reassembly",
    "annotation",
    "scaffolding",
    "mge",
]

MODULE_DIRS = {
    "preprocessing_sg": "1_preprocessing",
    "preprocessing_hic": "1_preprocessing",
    "assembly": "2_assembly",
    "alignment": "3_alignment",
    "coverage": "4_coverage",
    "contact": "5_contact",
    "binning": "6_binning",
    "reassembly": "7_reassembly",
    "annotation": "9_annotation",
    "scaffolding": "8_scaffolding",
    "mge": "10_MGE",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Nextflow driver for existing METAHI-T module commands."
    )
    parser.add_argument("--project-path", required=True)
    parser.add_argument("--sample-json", required=True)
    parser.add_argument("--modules", default="all")
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--result-name", default="results")
    parser.add_argument("--threads", type=int, default=80)
    parser.add_argument("--default-params", default="")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument(
        "--chain",
        action="store_true",
        help="Prefer outputs from earlier modules when their expected files exist.",
    )
    return parser.parse_args()


def parse_scalar(value):
    value = value.strip()
    if value in {"", "''", '""'}:
        return ""
    if value.lower() in {"null", "none", "~"}:
        return None
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if (
        (value.startswith('"') and value.endswith('"'))
        or (value.startswith("'") and value.endswith("'"))
    ):
        return value[1:-1]
    try:
        if any(char in value for char in (".", "e", "E")):
            return float(value)
        return int(value)
    except ValueError:
        return value


def parse_simple_yaml(path):
    root = {}
    stack = [(-1, root)]
    with open(path) as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.lstrip().startswith("#"):
                continue
            indent = len(raw_line) - len(raw_line.lstrip(" "))
            line = raw_line.split("#", 1)[0].rstrip()
            if ":" not in line:
                continue
            key, value = line.strip().split(":", 1)
            while stack and indent <= stack[-1][0]:
                stack.pop()
            current = stack[-1][1]
            if value.strip() == "":
                child = {}
                current[key] = child
                stack.append((indent, child))
            else:
                current[key] = parse_scalar(value)
    return root


def load_default_params(path, project_path=None):
    if not path:
        return {}
    config_path = Path(path)
    if not config_path.is_absolute() and not config_path.exists() and project_path is not None:
        candidate = Path(project_path) / config_path
        if candidate.exists():
            config_path = candidate
    if not config_path.exists():
        raise SystemExit(f"[ERROR] Default parameter file not found: {path}")
    try:
        import yaml  # type: ignore

        with open(config_path) as handle:
            data = yaml.safe_load(handle) or {}
    except ImportError:
        data = parse_simple_yaml(config_path)
    if not isinstance(data, dict):
        raise SystemExit(f"[ERROR] Default parameter file must contain a mapping: {path}")
    return data


def load_sample(path):
    with open(path) as handle:
        row = json.load(handle)
    row = {str(k).strip(): "" if v is None else str(v).strip() for k, v in row.items()}
    sample = row.get("sample", "").strip()
    if not sample:
        raise SystemExit("[ERROR] Samplesheet row is missing a non-empty sample value.")
    return row


def resolve_modules(text):
    text = (text or "all").strip()
    if text in {"all", "default", "short_read_all"}:
        return MODULE_ORDER
    modules = [item.strip() for item in text.split(",") if item.strip()]
    unknown = [item for item in modules if item not in MODULE_ORDER]
    if unknown:
        raise SystemExit(f"[ERROR] Unknown module name(s): {', '.join(unknown)}")
    return modules


def nonempty(row, *keys):
    for key in keys:
        value = row.get(key, "")
        if value:
            return value
    return ""


def as_path(value, project_path):
    if not value:
        return ""
    path = Path(value)
    if path.is_absolute():
        return str(path)
    return str((project_path / path).resolve())


def module_outdir(out_root, sample, module, result_name):
    if module == "preprocessing_sg":
        leaf = "results_sg" if result_name == "results" else f"{result_name}_sg"
    elif module == "preprocessing_hic":
        leaf = "results_hic" if result_name == "results" else f"{result_name}_hic"
    else:
        leaf = result_name
    return out_root / sample / MODULE_DIRS[module] / leaf


def module_aux_outdir(out_root, sample, module, result_name, suffix):
    leaf = f"{result_name}_{suffix}" if result_name else f"results_{suffix}"
    return out_root / sample / MODULE_DIRS[module] / leaf


def require_values(module, values):
    missing = [name for name, value in values.items() if not value]
    if missing:
        raise SystemExit(
            f"[ERROR] Module {module} is missing required samplesheet field(s): "
            + ", ".join(missing)
        )


def append_if(cmd, flag, value):
    if value != "":
        cmd.extend([flag, str(value)])


SPECIAL_FLAGS = {
    "checkm_db": "--checkm_db",
    "checkm2_db": "--checkm2_db",
    "gtdbtk_db": "--gtdbtk-db",
    "genomad_db": "--genomad-db",
}

NEGATIVE_FLAGS = {
    ("annotation", "skip_ani_screen"): "--no-skip-ani-screen",
    ("mge", "genomad_cleanup"): "--no-genomad-cleanup",
}

CONFIG_ONLY_KEYS = {
    "description",
    "threads",
    "method",
    "assembler",
}


def flag_for_key(key):
    return SPECIAL_FLAGS.get(key, f"--{key.replace('_', '-')}")


def module_defaults(default_params, module):
    modules = default_params.get("modules", {})
    if not isinstance(modules, dict):
        return {}
    section = modules.get(module, {})
    return section if isinstance(section, dict) else {}


def default_value(default_params, module, key):
    section = module_defaults(default_params, module)
    value = section.get(key, None)
    return "" if value is None else str(value)


def config_args(default_params, module):
    args = []
    section = module_defaults(default_params, module)
    for key, value in section.items():
        if key in CONFIG_ONLY_KEYS or value is None:
            continue
        if value == "":
            continue
        if key == "extra_args":
            if value:
                args.extend(shlex.split(str(value)))
            continue
        if isinstance(value, bool):
            if value:
                args.append(flag_for_key(key))
            else:
                negative = NEGATIVE_FLAGS.get((module, key))
                if negative:
                    args.append(negative)
            continue
        args.extend([flag_for_key(key), str(value)])
    return args


def extra_args(row, module, default_params=None):
    args = config_args(default_params or {}, module)
    text = nonempty(row, f"{module}_extra_args")
    if text:
        args.extend(shlex.split(text))
    return args


def preprocessing_project_path(project_path):
    direct = project_path / "external" / "bbmap" / "bbduk.sh"
    install = project_path / "installation" / "external" / "bbmap" / "bbduk.sh"
    if direct.exists():
        return project_path
    if install.exists():
        return project_path / "installation"
    return project_path


def subset_fastq_pair(read1, read2, outdir, prefix, read_count):
    outdir.mkdir(parents=True, exist_ok=True)
    out1 = outdir / f"{prefix}_R1.fastq.gz"
    out2 = outdir / f"{prefix}_R2.fastq.gz"
    if out1.exists() and out2.exists():
        return str(out1), str(out2)

    def copy_records(src, dest):
        with gzip.open(src, "rt") as inp, gzip.open(dest, "wt", compresslevel=6) as out:
            for line_no, line in enumerate(inp):
                if line_no >= read_count * 4:
                    break
                out.write(line)

    copy_records(read1, out1)
    copy_records(read2, out2)
    return str(out1), str(out2)


def stage_subset_inputs(row, project_path, out_root):
    read_count_text = nonempty(row, "preprocessing_subset_reads")
    if not read_count_text:
        return row

    read_count = int(read_count_text)
    if read_count <= 0:
        return row

    sample = row["sample"]
    stage_dir = out_root / sample / "_nextflow_inputs" / f"first_{read_count}_pairs"
    staged = dict(row)

    pairs = [
        ("sg1", "sg2", "preprocessing_sg1", "preprocessing_sg2", "shotgun"),
        ("hic1", "hic2", "preprocessing_hic1", "preprocessing_hic2", "hic"),
    ]
    for key1, key2, out_key1, out_key2, prefix in pairs:
        read1 = as_path(nonempty(row, key1), project_path)
        read2 = as_path(nonempty(row, key2), project_path)
        if read1 and read2:
            staged[out_key1], staged[out_key2] = subset_fastq_pair(
                read1, read2, stage_dir, prefix, read_count
            )
    return staged


def preprocessing_base(path):
    name = Path(path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if "_" in name:
        return name.rsplit("_", 1)[0]
    return name


def expected_preprocessed_pair(out_root, sample, result_name, source_path, module):
    outdir = module_outdir(out_root, sample, module, result_name)
    prefix = preprocessing_base(source_path)
    return (
        str(outdir / f"final_{prefix}_1.fastq.gz"),
        str(outdir / f"final_{prefix}_2.fastq.gz"),
    )


def use_existing_pair(pair, fallback_pair):
    if all(pair) and all(Path(item).exists() for item in pair):
        return pair
    return fallback_pair


def use_existing_file(path, fallback):
    if path and Path(path).exists():
        return str(path)
    return fallback


def use_chained_file(path, fallback):
    if path and (Path(path).exists() or not fallback):
        return str(path)
    return fallback


def chained_paths(module, row, project_path, out_root, result_name, sg1, sg2, hic1, hic2, fasta, bam, bin_dir):
    sample = row["sample"]
    use_preprocessed_reads = nonempty(row, "use_preprocessed_reads").lower() not in {
        "false", "0", "no", "n"
    }

    if use_preprocessed_reads and sg1 and sg2:
        sg1, sg2 = use_existing_pair(
            expected_preprocessed_pair(out_root, sample, result_name, sg1, "preprocessing_sg"),
            (sg1, sg2),
        )
    if use_preprocessed_reads and hic1 and hic2:
        hic1, hic2 = use_existing_pair(
            expected_preprocessed_pair(out_root, sample, result_name, hic1, "preprocessing_hic"),
            (hic1, hic2),
        )

    assembly_fasta = module_outdir(out_root, sample, "assembly", result_name) / "final_assembly.fasta"
    fasta = use_chained_file(assembly_fasta, fasta)

    alignment_bam = module_outdir(out_root, sample, "alignment", result_name) / "sorted_map.bam"
    bam = use_chained_file(alignment_bam, bam)

    binning_bins = module_outdir(out_root, sample, "binning", result_name) / "metahict" / "metahict_50_10_bins"
    bin_dir = use_chained_file(binning_bins, bin_dir)

    if module in {"annotation", "scaffolding"}:
        reassembly_bins = (
            module_outdir(out_root, sample, "reassembly", result_name)
            / "reassembled_bins"
        )
        bin_dir = use_chained_file(reassembly_bins, bin_dir)
        if module == "scaffolding" and reassembly_bins.exists():
            first_bin = sorted(reassembly_bins.glob("*.fa"))
            if first_bin:
                row["scaffold_fasta"] = str(first_bin[0])
        elif module == "scaffolding" and not row.get("scaffold_fasta", ""):
            row["scaffold_fasta"] = str(reassembly_bins / "bin1.fa")

    if module == "mge":
        reassembly_dir = module_outdir(out_root, sample, "reassembly", result_name)
        combined = reassembly_dir / "combined_contigs.fa"
        legacy_combined = reassembly_dir / "combined" / "combined_contigs.fa"
        if not combined.exists() and legacy_combined.exists():
            combined = legacy_combined
        if combined.exists() or not row.get("combined", ""):
            row["combined"] = str(combined)
        user_contact = row.get("contact", "")
        user_raw_contact = row.get("raw_contact", "")
        contact_dir = module_aux_outdir(out_root, sample, "contact", result_name, "MGE")
        normalized_contact = contact_dir / "denoised_contact_matrix_normcc.npz"
        raw_contact = contact_dir / "Raw_contact_matrix.npz"
        if not user_contact and not user_raw_contact:
            row["_prepare_mge_contact"] = "true"
            row["contact"] = str(normalized_contact)
            row["raw_contact"] = str(raw_contact)
        annotation_dir = module_outdir(out_root, sample, "annotation", result_name)
        host_taxonomy = annotation_dir / "classify" / "gtdbtk.bac120.summary.tsv"
        if host_taxonomy.exists() or not row.get("host_taxonomy", ""):
            row["host_taxonomy"] = str(host_taxonomy)

    return sg1, sg2, hic1, hic2, fasta, bam, bin_dir


def build_command(module, row, project_path, out_root, result_name, threads, chain=False, default_params=None):
    default_params = default_params or {}
    sample = row["sample"]
    metahict_py = str(project_path / "metahict.py")
    env_python = project_path / "conda_envs" / "metahict_env" / "bin" / "python"
    if os.environ.get("METAHICT_NF_CONTAINER") == "1":
        python = shutil.which("python") or sys.executable
    else:
        python = str(env_python) if env_python.exists() else (shutil.which("python") or sys.executable)
    outdir = str(module_outdir(out_root, sample, module, result_name))
    base = [python, metahict_py]

    sg1 = as_path(nonempty(row, "sg1", "shotgun_r1"), project_path)
    sg2 = as_path(nonempty(row, "sg2", "shotgun_r2"), project_path)
    hic1 = as_path(nonempty(row, "hic1", "hic_r1"), project_path)
    hic2 = as_path(nonempty(row, "hic2", "hic_r2"), project_path)
    fasta = as_path(nonempty(row, "fasta", "assembly", "reference"), project_path)
    bam = as_path(nonempty(row, "bam", "hic_bam"), project_path)
    bin_dir = as_path(nonempty(row, "bin_dir"), project_path)
    enzyme = nonempty(row, "enzyme") or "Sau3AI,MluCI"

    if chain:
        sg1, sg2, hic1, hic2, fasta, bam, bin_dir = chained_paths(
            module, row, project_path, out_root, result_name, sg1, sg2, hic1, hic2, fasta, bam, bin_dir
        )

    if module == "preprocessing_sg":
        pre_sg1 = as_path(nonempty(row, "preprocessing_sg1", "sg1", "shotgun_r1"), project_path)
        pre_sg2 = as_path(nonempty(row, "preprocessing_sg2", "sg2", "shotgun_r2"), project_path)
        require_values(module, {"sg1": pre_sg1, "sg2": pre_sg2})
        pre_project_path = preprocessing_project_path(project_path)
        return base + [
            "preprocessing", "-p", str(pre_project_path), "-1", pre_sg1, "-2", pre_sg2,
            "-o", outdir, "-t", str(threads),
        ] + extra_args(row, module, default_params)

    if module == "preprocessing_hic":
        pre_hic1 = as_path(nonempty(row, "preprocessing_hic1", "hic1", "hic_r1"), project_path)
        pre_hic2 = as_path(nonempty(row, "preprocessing_hic2", "hic2", "hic_r2"), project_path)
        require_values(module, {"hic1": pre_hic1, "hic2": pre_hic2})
        pre_project_path = preprocessing_project_path(project_path)
        return base + [
            "preprocessing", "-p", str(pre_project_path), "-1", pre_hic1, "-2", pre_hic2,
            "-o", outdir, "-t", str(threads), "--dedup",
        ] + extra_args(row, module, default_params)

    if module == "assembly":
        require_values(module, {"sg1": sg1, "sg2": sg2})
        cmd = base + [
            "assembly", "-p", str(project_path), "-1", sg1, "-2", sg2,
            "-o", outdir, "-t", str(threads),
        ]
        assembler = nonempty(row, "assembler") or default_value(default_params, module, "assembler")
        if assembler == "metaspades":
            cmd.append("--metaspades")
        elif assembler == "metaflye":
            cmd.append("--metaflye")
        else:
            cmd.append("--megahit")
        return cmd + extra_args(row, module, default_params)

    if module == "alignment":
        require_values(module, {"fasta": fasta, "hic1": hic1, "hic2": hic2})
        return base + [
            "alignment", "-p", str(project_path), "-r", fasta, "-1", hic1,
            "-2", hic2, "-o", outdir, "-t", str(threads),
        ] + extra_args(row, module, default_params)

    if module == "coverage":
        require_values(module, {"fasta": fasta, "sg1": sg1, "sg2": sg2})
        return base + [
            "coverage", "-p", str(project_path), "-r", fasta, "-1", sg1,
            "-2", sg2, "-o", outdir, "-t", str(threads),
        ] + extra_args(row, module, default_params)

    if module == "contact":
        require_values(module, {"fasta": fasta, "bam": bam})
        method = nonempty(row, "contact_method") or default_value(default_params, module, "method") or "normcc"
        return base + [
            "contact", method, "-p", str(project_path), "--bam", bam,
            "--fasta", fasta, "--out", outdir, "--enzyme", enzyme,
        ] + extra_args(row, module, default_params)

    if module == "binning":
        checkm_db = as_path(nonempty(row, "checkm_db"), project_path)
        require_values(module, {"fasta": fasta, "bam": bam})
        binning_project_path = project_path / "modules" if (project_path / "modules" / "6_binning" / "scripts").exists() else project_path
        cmd = base + [
            "binning", "--project_path", str(binning_project_path), "--fasta", fasta,
            "--bam", bam, "--output", outdir, "--enzyme", enzyme,
            "-t", str(threads), "--no-spades",
        ]
        append_if(cmd, "--checkm_db", checkm_db)
        return cmd + extra_args(row, module, default_params)

    if module == "reassembly":
        reassembly_project_path = project_path / "modules" if (project_path / "modules" / "7_reassembly" / "scripts").exists() else project_path
        require_values(
            module,
            {
                "bin_dir": bin_dir,
                "fasta": fasta,
                "hic1": hic1,
                "hic2": hic2,
                "sg1": sg1,
                "sg2": sg2,
                "bam": bam,
            },
        )
        return base + [
            "reassembly", "-p", str(reassembly_project_path), "--bin", bin_dir,
            "--assembly", fasta, "--hic1", hic1, "--hic2", hic2,
            "--sg1", sg1, "--sg2", sg2, "--bam", bam, "--outdir", outdir,
            "-t", str(threads),
        ] + extra_args(row, module, default_params)

    if module == "annotation":
        require_values(module, {"bin_dir": bin_dir})
        return base + [
            "annotation", "-p", str(project_path), "--bin", bin_dir,
            "--outdir", outdir, "-t", str(threads),
        ] + extra_args(row, module, default_params)

    if module == "scaffolding":
        scaffold_fasta = as_path(nonempty(row, "scaffold_fasta", "bin_fasta", "fasta"), project_path)
        require_values(module, {"scaffold_fasta": scaffold_fasta, "hic1": hic1, "hic2": hic2})
        cmd = base + [
            "scaffolding", "-p", str(project_path), "--fasta", scaffold_fasta,
            "--enzyme", enzyme, "--outdir", outdir, "--hic1", hic1,
            "--hic2", hic2, "-t", str(threads),
        ]
        if bam:
            cmd.extend(["--bam", bam])
        return cmd + extra_args(row, module, default_params)

    if module == "mge":
        combined = as_path(nonempty(row, "combined"), project_path)
        contact = as_path(nonempty(row, "contact"), project_path)
        raw_contact = as_path(nonempty(row, "raw_contact"), project_path)
        require_values(
            module,
            {
                "combined": combined,
                "contact": contact,
                "raw_contact": raw_contact,
            },
        )
        return base + [
            "mge", "-p", str(project_path), "--combined", combined,
            "--contact", contact, "--raw-contact", raw_contact,
            "--outdir", outdir,
            "-t", str(threads),
        ] + extra_args(row, module, default_params)

    raise AssertionError(f"Unhandled module: {module}")


def build_mge_prep_commands(row, project_path, out_root, result_name, threads, default_params=None):
    default_params = default_params or {}
    if row.get("_prepare_mge_contact", "").lower() != "true":
        return []

    sample = row["sample"]
    metahict_py = str(project_path / "metahict.py")
    env_python = project_path / "conda_envs" / "metahict_env" / "bin" / "python"
    if os.environ.get("METAHICT_NF_CONTAINER") == "1":
        python = shutil.which("python") or sys.executable
    else:
        python = str(env_python) if env_python.exists() else (shutil.which("python") or sys.executable)
    base = [python, metahict_py]

    combined = as_path(nonempty(row, "combined"), project_path)
    hic1 = as_path(nonempty(row, "hic1", "hic_r1"), project_path)
    hic2 = as_path(nonempty(row, "hic2", "hic_r2"), project_path)
    enzyme = nonempty(row, "enzyme") or "Sau3AI,MluCI"
    require_values("mge", {"combined": combined, "hic1": hic1, "hic2": hic2})

    alignment_outdir = str(module_aux_outdir(out_root, sample, "alignment", result_name, "MGE"))
    contact_outdir = str(module_aux_outdir(out_root, sample, "contact", result_name, "MGE"))
    mge_bam = str(Path(alignment_outdir) / "sorted_map.bam")

    return [
        base + [
            "alignment", "-p", str(project_path), "-r", combined, "-1", hic1,
            "-2", hic2, "-o", alignment_outdir, "-t", str(threads),
        ] + extra_args(row, "mge_alignment", default_params),
        base + [
            "contact", "normcc", "-p", str(project_path), "--bam", mge_bam,
            "--fasta", combined, "--out", contact_outdir, "--enzyme", enzyme,
        ] + extra_args(row, "mge_contact", default_params),
    ]


def run_command(cmd, dry_run):
    printable = " ".join(subprocess.list2cmdline([part]) for part in cmd)
    print(f"[METAHICT-NF] {printable}", flush=True)
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def main():
    args = parse_args()
    project_path = Path(args.project_path).resolve()
    out_root = Path(args.out_root).resolve()
    default_params = load_default_params(args.default_params, project_path)
    row = load_sample(args.sample_json)
    row = stage_subset_inputs(row, project_path, out_root)
    modules = resolve_modules(args.modules)

    for module in modules:
        outdir = module_outdir(out_root, row["sample"], module, args.result_name)
        if args.clean and outdir.exists():
            print(f"[METAHICT-NF] Removing previous output: {outdir}", flush=True)
            shutil.rmtree(outdir)
        if module == "mge" and args.clean:
            for prep_outdir in (
                module_aux_outdir(out_root, row["sample"], "alignment", args.result_name, "MGE"),
                module_aux_outdir(out_root, row["sample"], "contact", args.result_name, "MGE"),
            ):
                if prep_outdir.exists():
                    print(f"[METAHICT-NF] Removing previous output: {prep_outdir}", flush=True)
                    shutil.rmtree(prep_outdir)
        outdir.parent.mkdir(parents=True, exist_ok=True)
        cmd = build_command(
            module, row, project_path, out_root, args.result_name,
            args.threads, args.chain, default_params,
        )
        for prep_cmd in build_mge_prep_commands(
            row, project_path, out_root, args.result_name, args.threads, default_params
        ):
            run_command(prep_cmd, args.dry_run)
        run_command(cmd, args.dry_run)


if __name__ == "__main__":
    main()
