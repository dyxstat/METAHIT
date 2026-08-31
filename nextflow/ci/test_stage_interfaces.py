#!/usr/bin/env python3
"""Unit tests for Python stage interfaces replacing large shell drivers."""

from __future__ import annotations

import importlib.util
import ast
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]


def load_module(name: str, relative: str):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_plot_module(name: str, relative: str):
    """Load plot data functions without requiring matplotlib in source-only CI."""
    matplotlib = types.ModuleType("matplotlib")
    matplotlib.use = lambda _: None
    pyplot = types.ModuleType("matplotlib.pyplot")
    ticker = types.ModuleType("matplotlib.ticker")
    ticker.PercentFormatter = object
    with mock.patch.dict(
        sys.modules,
        {
            "matplotlib": matplotlib,
            "matplotlib.pyplot": pyplot,
            "matplotlib.ticker": ticker,
        },
    ):
        return load_module(name, relative)


MGE = load_module("metahict_mge_stage", "modules/10_MGE/mge.py")
PREPROCESSING = load_module(
    "metahict_preprocessing_stage",
    "modules/1_preprocessing/preprocessing.py",
)
ASSEMBLY = load_module("metahict_assembly_stage", "modules/2_assembly/assembly.py")
ALIGNMENT = load_module("metahict_alignment_stage", "modules/3_alignment/alignment.py")
COVERAGE = load_module("metahict_coverage_stage", "modules/4_coverage/coverage.py")
CONTACT = load_module("metahict_contact_stage", "modules/5_contact/contact.py")
ANNOTATION = load_module("metahict_annotation_stage", "modules/9_annotation/annotation.py")
REFINEMENT = load_module(
    "metahict_bin_refinement_stage",
    "modules/6_binning/bin_refinement.py",
)
BINNING_STAGE = load_module(
    "metahict_binning_stage",
    "modules/6_binning/run_binning.py",
)
REASSEMBLY = load_module(
    "metahict_reassemble_bins_stage",
    "modules/7_reassembly/reassemble_bins.py",
)
REASSEMBLY_DRIVER = load_module(
    "metahict_reassembly_driver",
    "modules/7_reassembly/run_reassembly.py",
)
BINNING_PLOT = load_plot_module(
    "metahict_binning_plot", "modules/6_binning/plot_binning_results.py"
)
REASSEMBLY_PLOT = load_plot_module(
    "metahict_reassembly_plot", "modules/7_reassembly/plot_reassembly.py"
)


class PreprocessingStageTest(unittest.TestCase):
    def test_selected_modules_use_one_preprocessing_stage_root(self) -> None:
        workflow = (ROOT / "nextflow/main_dsl2.nf").read_text()
        self.assertIn("params.preprocessing_dir = null", workflow)
        self.assertIn(
            'return directStageChannel("${params.preprocessing_dir}/${library}")',
            workflow,
        )
        self.assertNotIn("params.sg_preprocessing_dir", workflow)
        self.assertNotIn("params.hic_preprocessing_dir", workflow)

    def test_long_read_assembly_and_coverage_do_not_require_preprocessing(self) -> None:
        workflow = (ROOT / "nextflow/main_dsl2.nf").read_text()
        block = workflow.split("def shortReadShotgunStageChannel() {", 1)[1].split(
            "\ndef shotgunStageChannel()", 1
        )[0]
        self.assertLess(
            block.index(".filter { sample, row -> !rowLongReadType(row) }"),
            block.index("if (!params.preprocessing_dir)"),
        )

    def test_local_profile_caps_requests_to_detected_host_resources(self) -> None:
        configuration = (ROOT / "nextflow/nextflow.config").read_text()
        self.assertIn("process.resourceLimits", configuration)
        self.assertIn("cpus: params.local_resource_cpus", configuration)
        self.assertIn("memory: params.local_resource_memory", configuration)

    def test_nextflow_uses_yaml_module_maps_not_samplesheet_commands(self) -> None:
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        self.assertIn("def moduleArguments", processes)
        self.assertIn("def restrictionEnzyme", processes)
        self.assertNotRegex(processes, r"row\.[a-z_]+_extra_args")

    def test_prefix_handles_common_fastq_suffixes(self) -> None:
        self.assertEqual(PREPROCESSING.output_prefix(Path("sample_R1.fastq.gz"), None), "sample")
        self.assertEqual(PREPROCESSING.output_prefix(Path("sample_1.fq"), "library"), "library")

    def test_parser_preserves_documented_defaults(self) -> None:
        args = PREPROCESSING.build_parser().parse_args(
            ["-p", "/project", "-1", "r1.fq", "-2", "r2.fq", "-o", "out"]
        )
        self.assertEqual(args.threads, 8)
        self.assertEqual((args.minlen, args.trimq, args.k), (50, 10, 23))
        self.assertFalse(args.dedup)
        self.assertEqual(PREPROCESSING.DEFAULT_XMX, "25g")

        assembly_args = ASSEMBLY.build_parser().parse_args(
            ["-p", "/project", "-1", "r1.fq", "-2", "r2.fq", "-o", "out"]
        )
        alignment_args = ALIGNMENT.build_parser().parse_args(
            ["-p", "/project", "-r", "assembly.fa", "-1", "r1.fq", "-2", "r2.fq", "-o", "out"]
        )
        coverage_args = COVERAGE.build_parser().parse_args(
            ["-p", "/project", "-r", "assembly.fa", "-1", "r1.fq", "-2", "r2.fq", "-o", "out"]
        )
        reassembly_args = REASSEMBLY.build_parser().parse_args(
            ["-b", "bins", "-o", "out", "-1", "r1.fq", "-2", "r2.fq"]
        )
        annotation_args = ANNOTATION.build_parser().parse_args(
            ["-p", "/project", "--mag-dir", "mags", "--outdir", "out"]
        )
        mge_args = MGE.build_parser().parse_args(
            [
                "-p", "/project", "--fasta", "assembly.fa", "--host-dir", "hosts",
                "--contact", "contact.npz", "--raw-contact", "raw.npz", "--outdir", "out",
            ]
        )
        self.assertEqual((assembly_args.threads, assembly_args.memory), (16, 51))
        self.assertEqual(alignment_args.threads, 16)
        self.assertEqual((coverage_args.threads, coverage_args.memory), (16, "25g"))
        self.assertEqual((reassembly_args.threads, reassembly_args.memory), (16, 51))
        self.assertEqual(annotation_args.threads, 8)
        self.assertEqual(mge_args.threads, 16)

        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        block = processes.split("process PREPROCESSING {", 1)[1].split("\nprocess ", 1)[0]
        self.assertIn(
            "memory { moduleMemory('preprocessing', '32 GB') }",
            block,
        )
        self.assertIn("cpus { moduleThreads('preprocessing', 8) }", block)
        self.assertNotIn("time { moduleResource", processes)
        self.assertIn("modulePathArguments(\n        'preprocessing', ['trimming']", block)
        self.assertIn("['libraries', libraryKey]", block)
        self.assertIn("'deduplicate'", block)
        self.assertIn("library == 'hic'", block)
        self.assertIn("task.memory.toGiga() * 0.8", block)
        self.assertIn('--xmx "${memoryGb}g"', block)
        source = (ROOT / "modules/1_preprocessing/preprocessing.py").read_text()
        self.assertIn("run([clumpify, *common", source)

        analytical_defaults = [
            line
            for line in processes.splitlines()
            if "memory { moduleMemory" in line
        ]
        self.assertEqual(len(analytical_defaults), 12)
        self.assertIn("params.memory as nextflow.util.MemoryUnit", processes)
        self.assertIn("return override", processes)
        self.assertNotIn("Math.min(configured", processes)
        self.assertNotIn("params.memory_cap", processes)
        self.assertNotRegex(processes, r"\b(?:for|while)\s*\(")

        configuration = (ROOT / "nextflow/assets/metahict_configuration.yaml").read_text()
        resources = configuration.split("resources:\n", 1)[1].split("\nmodules:\n", 1)[0]
        expected_resources = {
            "preprocessing": (8, "32 GB"),
            "assembly": (16, "64 GB"),
            "alignment": (16, "32 GB"),
            "coverage": (16, "32 GB"),
            "contact": (1, "32 GB"),
            "binning": (16, "64 GB"),
            "reassembly": (16, "64 GB"),
            "scaffolding": (8, "32 GB"),
            "annotation": (8, "64 GB"),
            "mge": (16, "32 GB"),
        }
        for module, (threads, memory) in expected_resources.items():
            self.assertIn(
                f'  {module}: {{threads: {threads}, memory: "{memory}"}}',
                resources,
            )
        self.assertNotIn("mge_alignment:", resources)
        self.assertNotIn("mge_contact:", resources)
        self.assertIn("process MGE_CONTACT {", processes)
        self.assertIn("\n    cpus 1\n", processes)

    def test_dedup_switches_use_the_last_explicit_choice(self) -> None:
        required = ["-p", "/project", "-1", "r1.fq", "-2", "r2.fq", "-o", "out"]
        enabled = PREPROCESSING.build_parser().parse_args([*required, "--dedup"])
        disabled = PREPROCESSING.build_parser().parse_args(
            [*required, "--dedup", "--no-dedup"]
        )

        self.assertTrue(enabled.dedup)
        self.assertFalse(disabled.dedup)

    def test_adapter_reference_uses_conda_share_not_resolved_opt_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "metahict_env"
            executable_dir = prefix / "bin"
            implementation_dir = prefix / "opt" / "bbmap-39.10-0"
            adapter = prefix / "share" / "bbmap" / "resources" / "adapters.fa"
            executable_dir.mkdir(parents=True)
            implementation_dir.mkdir(parents=True)
            adapter.parent.mkdir(parents=True)
            adapter.write_text(">adapter\nACGT\n")
            implementation = implementation_dir / "bbdukOld.sh"
            implementation.touch()
            linked_executable = executable_dir / "bbdukOld.sh"
            linked_executable.symlink_to(implementation)

            observed = PREPROCESSING.bbtools_adapter_reference(
                str(linked_executable), environment={}
            )

            self.assertEqual(observed.resolve(), adapter.resolve())

    def test_adapter_reference_supports_upstream_bbtools_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "bbmap"
            executable = root / "bbduk.sh"
            adapter = root / "resources" / "adapters.fa"
            adapter.parent.mkdir(parents=True)
            executable.touch()
            adapter.write_text(">adapter\nACGT\n")

            observed = PREPROCESSING.bbtools_adapter_reference(
                str(executable), environment={}
            )

        self.assertEqual(observed.resolve(), adapter.resolve())


class BinningOutputLayoutTest(unittest.TestCase):
    @staticmethod
    def make_refinement_work(output: Path) -> Path:
        work = output / "work_files"
        final_source = work / "binsO"
        final_source.mkdir(parents=True)
        (final_source / "bin.1.fa").write_text(">contig1\nACGT\n")
        (work / "binsO.stats").write_text(
            "bin_id\tcompleteness\tcontamination\n"
            "bin.1\t90\t2\n"
        )
        candidate = work / "binsA"
        candidate.mkdir()
        (candidate / "candidate.fa").write_text(">candidate\nACGT\n")
        return work

    def test_default_layout_removes_refinement_internals(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            work = self.make_refinement_work(output)

            final_bins = REFINEMENT.publish_final_outputs(
                output, work, "binsO", keep_intermediates=False
            )

            self.assertEqual(final_bins, output / "final_bins")
            self.assertTrue((final_bins / "bin.1.fa").is_file())
            self.assertTrue((output / "final_bins_quality.tsv").is_file())
            self.assertFalse((output / "work_files").exists())
            self.assertFalse((output / "intermediates").exists())
            for obsolete in ("BIN", "fasta", "FINAL_BIN", "metahict_50_10_bins"):
                self.assertFalse((output / obsolete).exists())

    def test_requested_intermediates_have_one_parent_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            work = self.make_refinement_work(output)

            REFINEMENT.publish_final_outputs(
                output, work, "binsO", keep_intermediates=True
            )

            retained = output / "intermediates" / "refinement"
            self.assertTrue((retained / "binsO" / "bin.1.fa").is_file())
            self.assertTrue((retained / "binsA" / "candidate.fa").is_file())
            self.assertFalse((output / "work_files").exists())

    def test_run_parameter_manifest_uses_only_standard_python_types(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "run_parameters.yaml"
            args = types.SimpleNamespace(
                threads=32,
                enzyme="Sau3AI,MluCI",
                keep_temp=False,
                seed=None,
            )

            BINNING_STAGE.write_run_parameters(manifest, args)

            text = manifest.read_text()
            self.assertIn('metahict_version: "1.1.0"', text)
            self.assertIn("  threads: 32", text)
            self.assertIn('  enzyme: "Sau3AI,MluCI"', text)
            self.assertIn("  keep_temp: false", text)
            self.assertIn("  seed: null", text)


class PublicationContractTest(unittest.TestCase):
    def test_public_stage_parsers_explain_every_argument(self) -> None:
        stage_sources = (
            "modules/1_preprocessing/preprocessing.py",
            "modules/2_assembly/assembly.py",
            "modules/3_alignment/alignment.py",
            "modules/4_coverage/coverage.py",
            "modules/5_contact/contact.py",
            "modules/6_binning/run_binning.py",
            "modules/7_reassembly/run_reassembly.py",
            "modules/8_scaffolding/scaffolding.py",
            "modules/9_annotation/annotation.py",
            "modules/10_MGE/mge.py",
        )
        for relative in stage_sources:
            with self.subTest(stage=relative):
                tree = ast.parse((ROOT / relative).read_text())
                arguments = [
                    node
                    for node in ast.walk(tree)
                    if isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "add_argument"
                ]
                self.assertTrue(arguments)
                for argument in arguments:
                    self.assertTrue(
                        any(keyword.arg == "help" for keyword in argument.keywords),
                        f"{relative}:{argument.lineno} has no help text",
                    )

    def test_resolved_resources_are_forwarded_and_recorded(self) -> None:
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()

        def block(name: str) -> str:
            return processes.split(f"process {name} {{", 1)[1].split("\nprocess ", 1)[0]

        for name in (
            "PREPROCESSING", "ASSEMBLY", "ALIGNMENT", "COVERAGE", "BINNING",
            "REASSEMBLY", "ANNOTATION", "SCAFFOLDING", "MGE_ALIGNMENT", "MGE",
        ):
            self.assertIn('-t "${task.cpus}"', block(name), name)
        self.assertIn("cpus { moduleThreads('contact', 1) }", block("CONTACT"))
        self.assertIn("cpus 1", block("MGE_CONTACT"))

        memory_flags = {
            "PREPROCESSING": '--xmx "${memoryGb}g"',
            "ASSEMBLY": '--memory "${memoryGb}"',
            "COVERAGE": '--memory "${memoryGb}g"',
            "REASSEMBLY": '-m "${memoryGb}"',
        }
        for name, expected in memory_flags.items():
            self.assertIn("task.memory.toGiga() * 0.8", block(name))
            self.assertIn(expected, block(name))

        script_resource_use = {
            "modules/1_preprocessing/preprocessing.py": (
                'f"threads={args.threads}"', 'common = [f"-Xmx{xmx}"]'
            ),
            "modules/2_assembly/assembly.py": (
                '"-t", str(args.threads)', '"-m", str(memory)'
            ),
            "modules/3_alignment/alignment.py": (
                '"-t", str(args.threads)', '"-m", args.sort_memory'
            ),
            "modules/4_coverage/coverage.py": (
                'f"threads={args.threads}"', 'f"-Xmx{args.memory}"'
            ),
            "modules/6_binning/run_binning.py": (
                '"--threads", args.threads',
            ),
            "modules/7_reassembly/reassemble_bins.py": (
                "args.threads, args.memory",
            ),
            "modules/8_scaffolding/scaffolding.py": (
                '"-@", str(args.t)',
            ),
            "modules/9_annotation/annotation.py": (
                '"--cpus", str(args.threads)',
            ),
            "modules/10_MGE/mge.py": (
                '"--threads", args.threads', '"--ncpus", args.threads'
            ),
        }
        for relative, expected_fragments in script_resource_use.items():
            source = (ROOT / relative).read_text()
            for expected in expected_fragments:
                self.assertIn(expected, source, relative)

        self.assertEqual(processes.count("[INFO] Resources: threads="), 12)
        self.assertNotIn("resources.txt", processes)
        trace_config = (ROOT / "nextflow/nextflow.config").read_text()
        self.assertIn("exit,cpus,memory,submit", trace_config)

    def test_published_stage_directories_are_flat(self) -> None:
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        expected = (ROOT / "nextflow/tests/expected/workflow_stub_outputs.tsv").read_text()

        self.assertIn('${params.out_root}/${sample}/1_preprocessing"', processes)
        self.assertIn("filename == 'preprocessing' ? library", processes)
        self.assertIn("1_preprocessing/sg/final_sg_1.fastq.gz", expected)
        self.assertNotIn("1_preprocessing/sg/preprocessing/", expected)

        mappings = {
            "assembly": "2_assembly",
            "alignment": "3_alignment",
            "coverage": "4_coverage",
            "contact": "5_contact",
            "binning": "6_binning",
            "reassembly": "7_reassembly",
            "annotation": "9_annotation",
        }
        for internal, published in mappings.items():
            self.assertIn(
                f"filename == '{internal}' ? '{published}'",
                processes,
            )
        self.assertIn(
            "filename == 'scaffolding' ? \"8_scaffolding/${bin_id}\"",
            processes,
        )
        self.assertEqual(
            processes.count(
                "filename == 'mge' ? '10_MGE'"
            ),
            1,
        )
        self.assertEqual(processes.count("path('mge'), emit: results"), 1)
        self.assertNotIn("path('mge/*'), emit: results", processes)
        self.assertIn("10_MGE/intermediates/${stage}", processes)
        self.assertIn("return null", processes)

        duplicated_paths = (
            "2_assembly/assembly/",
            "3_alignment/alignment/",
            "4_coverage/coverage/",
            "5_contact/contact/",
            "6_binning/binning/",
            "7_reassembly/reassembly/",
            "8_scaffolding/scaffolding/",
            "9_annotation/annotation/",
            "10_MGE/mge_alignment/",
            "10_MGE/mge_contact/",
            "10_MGE/mge/",
        )
        for duplicated in duplicated_paths:
            self.assertNotIn(duplicated, expected)

    def test_explicit_stage_inputs_accept_current_and_legacy_layouts(self) -> None:
        workflow = (ROOT / "nextflow/main_dsl2.nf").read_text()

        self.assertEqual(workflow.count("def resolveStageDirectory"), 1)
        self.assertIn("markerSpec instanceof Collection ? markerSpec : [markerSpec]", workflow)
        self.assertNotIn("module_input_root", workflow)
        self.assertNotIn("def stageChannel", workflow)
        self.assertIn("directStageChannel(params.assembly_dir, 'final_assembly.fasta', ['assembly'])", workflow)
        self.assertIn("directStageChannel(params.alignment_dir, 'sorted_map.bam', ['alignment'])", workflow)
        self.assertIn(
            "directStageChannel(params.binning_dir, 'metahict/final_bins', ['binning'])",
            workflow,
        )
        self.assertIn("['alignment', 'mge_alignment', 'mge_alignment/mge_alignment']", workflow)
        self.assertIn("['contact', 'mge_contact', 'mge_contact/mge_contact']", workflow)

    def test_checkm2_processes_expose_companion_executables(self) -> None:
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        for process_name in ("BINNING", "REASSEMBLY", "SCAFFOLDING"):
            block = processes.split(f"process {process_name} {{", 1)[1].split("\nprocess ", 1)[0]
            self.assertIn("${params.conda_envs_path}/checkm2/bin", block)

    def test_scaffolding_uses_the_explicit_bin_and_realigns_hic(self) -> None:
        workflow = (ROOT / "nextflow/main_dsl2.nf").read_text()
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        block = processes.split("process SCAFFOLDING {", 1)[1].split("\nprocess ", 1)[0]

        self.assertIn("params.scaffolding_bin", workflow)
        self.assertIn("scaffoldingBinStageChannel()", workflow)
        self.assertIn("path(bin_fasta)", block)
        self.assertIn('--fasta "${bin_fasta}"', block)
        self.assertNotIn("reassembly_dir", block)
        self.assertNotIn("alignment_dir", block)
        self.assertIn('def bamArg = use_bam ? "--bam', block)
        self.assertIn("val(use_bam)", block)
        self.assertIn("params.scaffolding_bam", workflow)
        self.assertTrue((ROOT / "nextflow/assets/NO_SCAFFOLDING_BAM").is_file())
        self.assertNotIn("row.enzyme ?: 'Sau3AI,MluCI'", processes)

    def test_scaffolding_contact_path_is_metacc_only(self) -> None:
        source = (ROOT / "modules/8_scaffolding/scaffolding.py").read_text()
        heatmap = (ROOT / "modules/8_scaffolding/heatmap.py").read_text()

        self.assertFalse((ROOT / "modules/8_scaffolding/contact_matrix.py").exists())
        self.assertIn('os.path.join(args.p, "5_contact", "raw_contact.py")', source)
        self.assertNotIn("bin3c", source.lower())
        self.assertIn('"--contig-info", contig_csv', source)
        self.assertNotIn("obj.order", heatmap)
        self.assertIn('contig_info["name"]', heatmap)
        self.assertIn("Path(__file__).resolve().parents[2]", heatmap)

    def test_scaffolding_alignment_options_use_argparse_safe_assignments(self) -> None:
        source = (ROOT / "modules/8_scaffolding/scaffolding.py").read_text()
        self.assertIn('f"--bwa-options={args.bwa_options}"', source)
        self.assertIn('f"--samtools-filter={args.samtools_filter}"', source)
        self.assertNotIn('"--bwa-options", args.bwa_options', source)
        self.assertNotIn('"--samtools-filter", args.samtools_filter', source)

        parsed = ALIGNMENT.build_parser().parse_args(
            [
                "-p", "/project",
                "-r", "reference.fa",
                "-1", "hic_1.fastq.gz",
                "-2", "hic_2.fastq.gz",
                "-o", "alignment",
                "--bwa-options=-5SP",
                "--samtools-filter=-F 0x900",
            ]
        )
        self.assertEqual(parsed.bwa_options, "-5SP")
        self.assertEqual(parsed.samtools_filter, "-F 0x900")

    def test_scaffolding_validates_optional_bam_reference(self) -> None:
        source = (ROOT / "modules/8_scaffolding/scaffolding.py").read_text()
        self.assertIn("def validate_bam_reference", source)
        self.assertIn("validate_bam_reference(args.bam, filt_fa)", source)

    def test_scaffolding_separates_results_from_intermediates(self) -> None:
        source = (ROOT / "modules/8_scaffolding/scaffolding.py").read_text()
        real_expected = (
            ROOT / "nextflow/tests/expected/example_dataset_outputs.tsv"
        ).read_text()
        stub_expected = (
            ROOT / "nextflow/tests/expected/workflow_stub_outputs.tsv"
        ).read_text()

        self.assertIn('os.path.join(args.outdir, "scaffolded_bin.fa")', source)
        self.assertIn('os.path.join(args.outdir, "scaffolded_bin.agp")', source)
        self.assertIn('os.path.join(args.outdir, "intermediates")', source)
        self.assertIn("cleanup_intermediates(intermediate_dir)", source)
        self.assertIn('os.path.join(args.outdir, "quality")', source)
        self.assertIn("8_scaffolding/*/scaffolded_bin.fa", real_expected)
        self.assertIn("8_scaffolding/bin1/scaffolded_bin.fa", stub_expected)
        self.assertNotIn("8_scaffolding/*/scaffolds.fa", real_expected)

    def test_every_scaffolding_parser_argument_is_used(self) -> None:
        source = (ROOT / "modules/8_scaffolding/scaffolding.py").read_text()
        tree = ast.parse(source)
        declared = set()
        used = set()

        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"
            ):
                flags = [
                    arg.value
                    for arg in node.args
                    if isinstance(arg, ast.Constant) and isinstance(arg.value, str)
                ]
                destination = next(
                    (
                        keyword.value.value
                        for keyword in node.keywords
                        if keyword.arg == "dest"
                        and isinstance(keyword.value, ast.Constant)
                    ),
                    None,
                )
                if destination is None:
                    preferred = next((flag for flag in flags if flag.startswith("--")), flags[0])
                    destination = preferred.lstrip("-").replace("-", "_")
                declared.add(destination)
            if (
                isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id == "args"
            ):
                used.add(node.attr)
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "getattr"
                and len(node.args) >= 2
                and isinstance(node.args[0], ast.Name)
                and node.args[0].id == "args"
                and isinstance(node.args[1], ast.Constant)
            ):
                used.add(node.args[1].value)

        self.assertEqual(declared - used, set())
        self.assertIn('"--resolution"', source)
        self.assertNotIn('parser.add_argument("-m"', source)

    def test_reassembly_uses_nextflow_memory_allocation(self) -> None:
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        block = processes.split("process REASSEMBLY {", 1)[1].split("\nprocess ", 1)[0]
        self.assertIn("task.memory.toGiga() * 0.8", block)
        self.assertIn('-m "${memoryGb}"', block)

    def test_seqtk_output_is_compressed_through_python(self) -> None:
        source = (ROOT / "modules/7_reassembly/select_reassembly_reads.py").read_text()
        self.assertIn("shutil.copyfileobj(process.stdout, output", source)
        self.assertNotIn("stdout=output", source)


class AssemblyStageTest(unittest.TestCase):
    def test_parser_defaults_to_megahit(self) -> None:
        args = ASSEMBLY.build_parser().parse_args(
            ["-p", "/project", "-1", "r1.fq", "-2", "r2.fq", "-o", "out"]
        )
        self.assertEqual(args.assembler, "megahit")

    def test_parser_accepts_one_typed_long_read_file(self) -> None:
        args = ASSEMBLY.build_parser().parse_args(
            [
                "-p", "/project", "--long-reads", "reads.fastq.gz",
                "--long-read-type", "pacbio-hifi", "--metaflye", "-o", "out",
            ]
        )
        self.assertEqual(args.long_read_type, "pacbio-hifi")
        self.assertEqual(args.assembler, "metaflye")
        self.assertIsNone(args.r2)

    def test_invalid_kmer_configuration_is_rejected(self) -> None:
        args = ASSEMBLY.build_parser().parse_args(
            ["-p", "/project", "-1", "r1.fq", "-2", "r2.fq", "-o", "out", "--k-step", "11"]
        )
        with self.assertRaises(ValueError):
            ASSEMBLY.validate_kmers(args)


class AlignmentStageTest(unittest.TestCase):
    def test_cigar_match_length_counts_aligned_bases(self) -> None:
        self.assertEqual(ALIGNMENT.match_length("10M2I5=3X4S"), 18)
        self.assertEqual(ALIGNMENT.match_length("*"), 0)


class CoverageStageTest(unittest.TestCase):
    def test_parser_accepts_outdir_alias_and_numeric_filters(self) -> None:
        args = COVERAGE.build_parser().parse_args(
            ["-p", "/project", "-1", "r1", "-2", "r2", "-r", "ref", "--outdir", "out",
             "--percent-identity", "98.5", "--min-mapq", "10"]
        )
        self.assertEqual((args.percent_identity, args.min_mapq), (98.5, 10))

    def test_parser_accepts_one_long_read_file(self) -> None:
        args = COVERAGE.build_parser().parse_args(
            ["-p", "/project", "--long-reads", "reads.fastq.gz", "-r", "ref", "-o", "out"]
        )
        self.assertEqual(args.long_reads, Path("reads.fastq.gz"))
        self.assertIsNone(args.reads2)


class ContactStageTest(unittest.TestCase):
    def test_parser_accepts_threshold_alias(self) -> None:
        args = CONTACT.build_parser().parse_args(
            ["normcc", "-p", "/project", "--bam", "map.bam", "--fasta", "assembly.fa",
             "--out", "out", "--enzyme", "MluCI", "--thres", "7.5"]
        )
        self.assertEqual(args.spurious_contact_percent, 7.5)


class AnnotationStageTest(unittest.TestCase):
    def test_help_uses_a_user_facing_scientific_description(self) -> None:
        help_text = ANNOTATION.build_parser().format_help()
        self.assertIn("Classify metagenome-assembled genomes (MAGs) with GTDB-Tk.", help_text)
        self.assertNotIn("shell orchestration", help_text)

    def test_command_uses_argument_vector_for_extra_options(self) -> None:
        args = ANNOTATION.build_parser().parse_args(
            ["-p", "/project", "--mag-dir", "mags", "--outdir", "out",
             "--gtdbtk-extra-args", "--genes --no_mash"]
        )
        command = ANNOTATION.build_command(args, "gtdbtk", Path("db"), Path("tmp"))
        self.assertEqual(command[-2:], ["--genes", "--no_mash"])

    def test_nextflow_annotation_consumes_the_explicit_mag_directory(self) -> None:
        workflow = (ROOT / "nextflow/main_dsl2.nf").read_text()
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        block = processes.split("process ANNOTATION {", 1)[1].split("\nprocess ", 1)[0]

        self.assertIn("params.annotation_mag_dir", workflow)
        self.assertIn("def annotationMagStageChannel()", workflow)
        self.assertIn("ANNOTATION(annotationMagStageChannel())", workflow)
        self.assertIn("path(mag_dir)", block)
        self.assertIn('--mag-dir "${mag_dir}"', block)
        self.assertNotIn("reassembly_dir", block)
        self.assertIn('file("${reassembly_directory}/reassembled_bins"', workflow)


class MgeStageTest(unittest.TestCase):
    def test_duplicate_fasta_ids_are_renamed_deterministically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.fa"
            destination = root / "unique.fa"
            mapping = root / "map.tsv"
            source.write_text(">contig description\nAAAA\n>contig\nCCCC\n>other\nGGGG\n")
            total, renamed = MGE.prepare_unique_fasta(source, destination, mapping)
            self.assertEqual((total, renamed), (3, 1))
            self.assertEqual(
                [line for line in destination.read_text().splitlines() if line.startswith(">")],
                [">contig description", ">contig__dup2", ">other"],
            )

    def test_parser_accepts_documented_genomad_alias(self) -> None:
        args = MGE.build_parser().parse_args(
            [
                "-p", "/project", "--fasta", "assembly.fa",
                "--host-dir", "hosts",
                "--contact", "contact.npz", "--raw-contact", "raw.npz",
                "--outdir", "out", "--genomad_db", "db",
            ]
        )
        self.assertEqual(args.genomad_db, "db")
        self.assertEqual(args.fasta, "assembly.fa")
        self.assertEqual(args.host_dir, "hosts")
        self.assertTrue(args.genomad_cleanup)

    def test_explicit_contig_to_host_map_supports_generic_and_metahict_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly = root / "assembly.fa"
            assembly.write_text(
                ">long_read_contig\nAAAA\n>bin9|reassembled_contig\nCCCC\n>unbinned\nGGGG\n"
            )
            hosts = root / "hosts"
            hosts.mkdir()
            (hosts / "long_read_mag.fa").write_text(">long_read_contig\nAAAA\n")
            (hosts / "bin9.fa").write_text(">reassembled_contig\nCCCC\n")
            mapping_file = root / "contig_to_host.tsv"

            mapping = MGE.build_contig_to_host(assembly, hosts, mapping_file)
            mapping_lines = mapping_file.read_text().splitlines()

        self.assertEqual(
            mapping,
            {
                "long_read_contig": "long_read_mag",
                "bin9|reassembled_contig": "bin9",
            },
        )
        self.assertEqual(
            mapping_lines,
            [
                "contig_id\thost_id",
                "long_read_contig\tlong_read_mag",
                "bin9|reassembled_contig\tbin9",
            ],
        )

    def test_mge_nextflow_inputs_are_not_reassembly_directories(self) -> None:
        workflow = (ROOT / "nextflow/main_dsl2.nf").read_text()
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        alignment = processes.split("process MGE_ALIGNMENT {", 1)[1].split("\nprocess ", 1)[0]
        contact = processes.split("process MGE_CONTACT {", 1)[1].split("\nprocess ", 1)[0]
        mge = processes.split("process MGE {", 1)[1].split("\nprocess ", 1)[0]

        for block in (alignment, contact, mge):
            self.assertNotIn("reassembly_dir", block)
        self.assertIn("path(fasta)", alignment)
        self.assertIn("path(fasta)", contact)
        self.assertIn("path(host_dir)", mge)
        self.assertIn('--fasta "${fasta}"', mge)
        self.assertIn('--host-dir "${host_dir}"', mge)
        self.assertNotIn("--mode", mge)
        self.assertNotIn("process MGE_DETECTION", processes)
        self.assertIn("/ccfind_env/bin", mge)
        self.assertIn("def mgeFastaStageChannel()", workflow)
        self.assertIn("def mgeHostStageChannel()", workflow)
        self.assertNotIn('contig.startswith("bin")', (ROOT / "modules/10_MGE/mge.py").read_text())

    def test_mge_ccfind_environment_exposes_ssearch36(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            binary = project / "conda_envs" / "ccfind_env" / "bin" / "ssearch36"
            binary.parent.mkdir(parents=True)
            binary.write_text("#!/bin/sh\nexit 0\n")
            binary.chmod(0o755)

            environment = MGE.ccfind_runtime_environment(project, {"PATH": ""})

        self.assertEqual(environment["PATH"].split(os.pathsep)[0], str(binary.parent))

    def test_mge_default_record_describes_one_complete_workflow(self) -> None:
        defaults = (ROOT / "nextflow/assets/metahict_configuration.yaml").read_text()
        mge_defaults = defaults.split("\n  mge:\n", 1)[1].split("\n  reassembly:\n", 1)[0]
        self.assertNotIn("mode:", mge_defaults)
        self.assertIn("keep_intermediates: false", mge_defaults)
        self.assertIn("\n    pairs:\n", mge_defaults)
        self.assertIn("\n    alignment:\n", mge_defaults)
        self.assertIn("\n    contact:\n", mge_defaults)

    def test_mge_publication_hides_supporting_files_by_default(self) -> None:
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        self.assertIn(
            "modulePathBoolean('mge', [], 'keep_intermediates', false)",
            processes,
        )
        self.assertIn("mgeIntermediatePublishPath(filename, 'alignment')", processes)
        self.assertIn("mgeIntermediatePublishPath(filename, 'contact')", processes)
        self.assertNotIn("10_MGE/alignment", processes)
        self.assertNotIn("10_MGE/contact", processes)
        self.assertNotIn("10_MGE/results", processes)

    def test_internal_mge_settings_come_from_one_module_map(self) -> None:
        processes = (ROOT / "nextflow/modules/local/metahict_modules.nf").read_text()
        self.assertIn("alignmentConfigurationArguments('mge', ['alignment'])", processes)
        self.assertIn("contactConfiguration('mge', ['contact'])", processes)
        self.assertIn("mgeAnalysisArguments()", processes)
        self.assertNotIn("moduleArguments('mge_alignment')", processes)
        self.assertNotIn("moduleSettings('mge_contact')", processes)


class RefinementStageTest(unittest.TestCase):
    def test_checkm2_environment_exposes_diamond(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            executable_dir = Path(temporary) / "checkm2" / "bin"
            executable_dir.mkdir(parents=True)
            checkm2 = executable_dir / "checkm2"
            diamond = executable_dir / "diamond"
            checkm2.touch(mode=0o755)
            diamond.touch(mode=0o755)
            with mock.patch.dict(os.environ, {"PATH": "/usr/bin"}):
                environment = REFINEMENT.checkm2_environment(str(checkm2))
        self.assertEqual(environment["PATH"].split(os.pathsep)[0], str(executable_dir))

    def test_input_bins_are_size_filtered_and_headers_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            destination = root / "binsA"
            source.mkdir()
            (source / "keep.fasta").write_text(">a=b\n" + "A" * 100 + "\n")
            (source / "small.fa").write_text(">small\nA\n")
            count = REFINEMENT.prepare_bin_set(source, destination, 10, 1000)
            self.assertEqual(count, 1)
            self.assertEqual((destination / "keep.fa").read_text().splitlines()[0], ">a_b")

    def test_quality_count_uses_documented_inclusive_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stats = Path(temporary) / "bins.stats"
            stats.write_text(
                "bin\tcompleteness\tcontamination\n"
                "a\t50\t10\n"
                "b\t49.9\t1\n"
                "c\t90\t10.1\n"
            )
            self.assertEqual(REFINEMENT.good_bin_count(stats, 50, 10), 1)

    def test_plot_scores_accept_decimal_thresholds_and_empty_results(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stats = Path(temporary) / "bins.stats"
            stats.write_text(
                "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
                "keep\t50.5\t10.25\t45\tUnknown\t1000\t50000\n"
                "incomplete\t50.4\t1\t45\tUnknown\t900\t50000\n"
                "contaminated\t90\t10.26\t45\tUnknown\t800\t50000\n"
            )
            completion, contamination = BINNING_PLOT.load_scores(stats, 50.5, 10.25)
            empty_completion, empty_contamination = BINNING_PLOT.load_scores(
                stats, 99.0, 0.0
            )
        self.assertEqual((completion, contamination), ([50.5], [10.25]))
        self.assertEqual((empty_completion, empty_contamination), ([], []))


class ReassemblyStageTest(unittest.TestCase):
    @staticmethod
    def make_public_results(output: Path, include_quality: bool = True) -> None:
        bins = output / "reassembled_bins"
        bins.mkdir(parents=True)
        (bins / "bin1.fa").write_text(">contig1\nACGT\n")
        (output / "combined_contigs.fa").write_text(">bin1|contig1\nACGT\n")
        (output / "residual_contigs.fa").write_text("")
        (output / "read_selection_summary.json").write_text(
            '{"selection": {"selected_pairs": 1}}\n'
        )
        if include_quality:
            (output / "reassembled_bins_quality.tsv").write_text(
                "bin\tcompleteness\tcontamination\n"
                "bin1\t90\t2\n"
            )

    def test_default_layout_removes_large_reassembly_intermediates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            self.make_public_results(output)
            (output / "hic.name_sorted.bam").write_bytes(b"bam")
            (output / "new_sg_forward.fastq.gz").write_bytes(b"reads")
            (output / "all_intra_insert_sizes.tsv.gz").write_bytes(b"table")
            (output / "input_assembly").mkdir()
            (output / "work_files").mkdir()
            (output / "reassembled_bins.checkm2").mkdir()
            (output / "reassembly_results.png").write_bytes(b"figure")
            args = types.SimpleNamespace(
                skip_checkm2=False,
                keep_temp=False,
                threads=32,
            )

            REASSEMBLY_DRIVER.finalize_output_layout(output, args)

            self.assertTrue((output / "reassembled_bins" / "bin1.fa").is_file())
            self.assertTrue((output / "combined_contigs.fa").is_file())
            self.assertTrue((output / "reassembled_bins_quality.tsv").is_file())
            self.assertTrue((output / "figures" / "reassembly_results.png").is_file())
            self.assertTrue((output / "run_parameters.yaml").is_file())
            for obsolete in (
                "hic.name_sorted.bam",
                "new_sg_forward.fastq.gz",
                "all_intra_insert_sizes.tsv.gz",
                "input_assembly",
                "work_files",
                "reassembled_bins.checkm2",
                "intermediates",
            ):
                self.assertFalse((output / obsolete).exists(), obsolete)

    def test_requested_reassembly_intermediates_are_grouped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            self.make_public_results(output, include_quality=False)
            (output / "readname_sg_in_hic.txt").write_text("read1\n")
            (output / "original_bins").mkdir()
            args = types.SimpleNamespace(
                skip_checkm2=True,
                keep_temp=True,
                threads=8,
            )

            REASSEMBLY_DRIVER.finalize_output_layout(output, args)

            self.assertTrue(
                (output / "intermediates" / "read_selection" / "readname_sg_in_hic.txt").is_file()
            )
            self.assertTrue(
                (output / "intermediates" / "reassembly" / "original_bins").is_dir()
            )
            self.assertTrue((output / "reassembled_bins" / "bin1.fa").is_file())

    def test_pipeline_process_arguments_are_normalized_to_strings(self) -> None:
        with mock.patch.object(REASSEMBLY.subprocess, "Popen") as popen:
            REASSEMBLY.start_process(
                ["filter", Path("bins"), 2, 5], stdin=REASSEMBLY.subprocess.PIPE
            )
        self.assertEqual(popen.call_args.args[0], ["filter", "bins", "2", "5"])

    def test_checkm2_environment_exposes_diamond(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            executable_dir = Path(temporary) / "checkm2" / "bin"
            executable_dir.mkdir(parents=True)
            checkm2 = executable_dir / "checkm2"
            diamond = executable_dir / "diamond"
            checkm2.touch(mode=0o755)
            diamond.touch(mode=0o755)
            with mock.patch.dict(os.environ, {"PATH": "/usr/bin"}):
                environment = REASSEMBLY.checkm2_environment(str(checkm2))
        self.assertEqual(environment["PATH"].split(os.pathsep)[0], str(executable_dir))

    def test_combined_contigs_have_bin_and_residual_prefixes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bins = root / "bins"
            bins.mkdir()
            (bins / "bin1.fa").write_text(">contig\nAAAA\n")
            residual = root / "residual.fa"
            residual.write_text(">contig\nCCCC\n")
            output = root / "combined.fa"
            REASSEMBLY.write_combined_contigs(bins, residual, output)
            headers = [line for line in output.read_text().splitlines() if line.startswith(">")]
            self.assertEqual(headers, [">bin1|contig", ">residual|contig"])

    def test_minimum_length_filter_is_sequence_based(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.fa"
            output = root / "output.fa"
            source.write_text(">short\nAAA\n>keep\nAAA\nAAA\n")
            count = REASSEMBLY.write_minimum_length_fasta(source, output, 5)
            self.assertEqual(count, 1)
            self.assertIn(">keep", output.read_text())

    def test_reassembly_plot_accepts_decimal_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stats = Path(temporary) / "reassembly.stats"
            stats.write_text(
                "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
                "keep\t50.5\t10.25\t45\tUnknown\t1000.0\t50000\n"
            )
            n50, completion, contamination = REASSEMBLY_PLOT.load_scores(
                stats, 50.5, 10.25
            )
        self.assertEqual((n50, completion, contamination), ([1000], [50.5], [10.25]))

    def test_parallel_spades_resources_respect_total_allocation(self) -> None:
        workers, cpus_per_job, memory_per_job = REASSEMBLY.parallel_spades_resources(
            100, 32, 102
        )
        self.assertEqual((workers, cpus_per_job, memory_per_job), (12, 2, 8))
        self.assertLessEqual(workers * cpus_per_job, 32)
        self.assertLessEqual(workers * memory_per_job, 102)

    def test_recruitment_resume_requires_completion_marker(self) -> None:
        source = (ROOT / "modules/7_reassembly/reassemble_bins.py").read_text()
        self.assertIn('reads_output / ".recruitment_complete"', source)
        self.assertIn('recruitment_marker.is_file()', source)


if __name__ == "__main__":
    unittest.main()
