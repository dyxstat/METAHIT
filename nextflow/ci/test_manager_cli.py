#!/usr/bin/env python3
"""Unit tests for the dependency-free METAHICT management CLI."""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import metahict_manager as manager


class LauncherTest(unittest.TestCase):
    def test_launcher_uses_conda_python(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            conda_base = root / "conda-base"
            fake_bin.mkdir()
            (conda_base / "bin").mkdir(parents=True)
            (conda_base / "bin" / "python").symlink_to(Path(sys.executable))
            fake_conda = fake_bin / "conda"
            fake_conda.write_text(
                "#!/bin/sh\n"
                f"printf '%s\\n' {shlex.quote(str(conda_base))}\n"
            )
            fake_conda.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"
            result = subprocess.run(
                [str(manager.PROJECT_ROOT / "metahict"), "--version"],
                check=True,
                capture_output=True,
                text=True,
                env=env,
            )
        self.assertEqual(result.stdout.strip(), f"METAHICT {manager.VERSION}")


class SourceMetadataTest(unittest.TestCase):
    def test_macos_archive_metadata_is_not_source(self) -> None:
        self.assertTrue(manager.is_macos_metadata_path(Path("modules/._stage.sh")))
        self.assertTrue(manager.is_macos_metadata_path(Path("__MACOSX/module.py")))
        self.assertTrue(manager.is_macos_metadata_path(Path("docs/.DS_Store")))
        self.assertFalse(manager.is_macos_metadata_path(Path("modules/stage.py")))


class ExplicitLockTest(unittest.TestCase):
    def test_headers_are_not_part_of_package_identity(self) -> None:
        lines = ["# platform: linux-64", "@EXPLICIT", "https://b", "", "https://a"]
        self.assertEqual(manager.explicit_urls(lines), ["https://a", "https://b"])

    def test_direct_python_dependencies_have_reportable_versions(self) -> None:
        wheel = (
            "Binning-refiner @ https://example.org/"
            "Binning_refiner-1.4.3-py3-none-any.whl#sha256=abc"
        )
        git = "bin3c-python3-metahict @ git+https://example.org/bin3c.git@181d80f"
        self.assertEqual(
            manager.pinned_python_requirement(wheel),
            ("Binning-refiner", "1.4.3", "wheel"),
        )
        self.assertEqual(
            manager.pinned_python_requirement(git),
            ("bin3c-python3-metahict", "181d80f", "git-commit"),
        )


class EnvironmentMigrationTest(unittest.TestCase):
    def test_old_nextflow_environment_is_used_as_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            environment_root = Path(temporary)
            legacy = environment_root / "metahict_venv"
            legacy.mkdir()
            with mock.patch.object(manager, "ENV_ROOT", environment_root):
                resolved = manager.environment_prefix("metahict_nextflow_env")
        self.assertEqual(resolved, legacy)

    def test_new_nextflow_environment_takes_precedence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            environment_root = Path(temporary)
            (environment_root / "metahict_venv").mkdir()
            current = environment_root / "metahict_nextflow_env"
            current.mkdir()
            with mock.patch.object(manager, "ENV_ROOT", environment_root):
                resolved = manager.environment_prefix("metahict_nextflow_env")
        self.assertEqual(resolved, current)


class NextflowCommandTest(unittest.TestCase):
    def setUp(self) -> None:
        self.runtime_patcher = mock.patch.object(manager, "verify_runtime")
        self.runtime_patcher.start()
        self.addCleanup(self.runtime_patcher.stop)

    def test_all_descriptive_entry_names_are_accepted(self) -> None:
        parser = manager.build_parser()
        for entry in manager.ENTRY_MODULES:
            with self.subTest(entry=entry):
                args = parser.parse_args(
                    [
                        "run",
                        "--entry-module",
                        entry,
                        "--samplesheet",
                        "samples.csv",
                        "--outdir",
                        "results",
                    ]
                )
                self.assertEqual(args.entry_module, entry)

    def test_command_has_documented_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample,sg1,sg2,hic1,hic2\n")
            args = argparse.Namespace(
                samplesheet=str(samplesheet),
                outdir=str(root / "results"),
                report_dir=None,
                work_dir=None,
                threads=12,
                memory="96 GB",
                checkm_db=None,
                checkm2_db=None,
                gtdbtk_db=None,
                genomad_db=None,
                resume=False,
            )
            command = manager.build_nextflow_command(
                args, resource_limits=(24, 120 * 1024 ** 3)
            )
        self.assertIn("nextflow/main_dsl2.nf", command[2])
        self.assertEqual(
            command[command.index("-params-file") + 1],
            str(manager.DEFAULT_CONFIG.resolve()),
        )
        self.assertEqual(command[command.index("-profile") + 1], "local")
        self.assertEqual(command[command.index("--entry_module") + 1], "all")
        self.assertEqual(command[command.index("--threads") + 1], "12")
        self.assertEqual(command[command.index("--memory") + 1], "96 GB")
        self.assertEqual(
            command[command.index("--local_resource_cpus") + 1], "24"
        )
        self.assertEqual(
            command[command.index("--local_resource_memory") + 1], "122880 MB"
        )
        self.assertIn("uniref100.KO.1.dmnd", command[command.index("--checkm2_db") + 1])

    def test_local_resource_warning_reports_automatic_caps(self) -> None:
        args = argparse.Namespace(
            config=str(manager.DEFAULT_CONFIG),
            entry_module="all",
            threads=None,
            memory=None,
        )
        warning = io.StringIO()
        with contextlib.redirect_stderr(warning):
            manager.warn_if_resources_are_capped(args, 8, 48 * 1024 ** 3)
        text = warning.getvalue()
        self.assertIn("up to 16 threads/64 GB requested", text)
        self.assertIn("8 threads/48 GB detected", text)
        self.assertIn("Nextflow will cap each task", text)

        no_warning = io.StringIO()
        with contextlib.redirect_stderr(no_warning):
            manager.warn_if_resources_are_capped(args, 64, 192 * 1024 ** 3)
        self.assertEqual(no_warning.getvalue(), "")

    def test_named_entry_and_upstream_inputs_are_forwarded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "alignment",
                    "--samplesheet",
                    str(root / "samples.csv"),
                    "--outdir",
                    str(root / "alignment-run"),
                    "--assembly-dir",
                    str(root / "assembly"),
                    "--preprocessing-dir",
                    str(root / "preprocessing"),
                ]
            )
            command = manager.build_nextflow_command(args)
        self.assertEqual(command[command.index("--entry_module") + 1], "alignment")
        self.assertNotIn("--checkm_db", command)
        self.assertNotIn("--checkm2_db", command)
        self.assertNotIn("--gtdbtk_db", command)
        self.assertNotIn("--genomad_db", command)
        self.assertEqual(
            command[command.index("--assembly_dir") + 1],
            str((root / "assembly").resolve()),
        )
        self.assertEqual(
            command[command.index("--preprocessing_dir") + 1],
            str((root / "preprocessing").resolve()),
        )

    def test_scaffolding_forwards_one_explicit_bin(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bin_fasta = root / "reassembled_bins" / "bin.1.fa"
            bin_fasta.parent.mkdir()
            bin_fasta.write_text(">contig\nACGT\n")
            bam = root / "bin.1.bam"
            bam.touch()
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "scaffolding",
                    "--samplesheet",
                    str(root / "samples.csv"),
                    "--outdir",
                    str(root / "scaffolding-run"),
                    "--scaffolding-bin",
                    str(bin_fasta),
                    "--scaffolding-bam",
                    str(bam),
                    "--preprocessing-dir",
                    str(root / "preprocessing"),
                ]
            )
            command = manager.build_nextflow_command(args)
        self.assertEqual(
            command[command.index("--scaffolding_bin") + 1],
            str(bin_fasta.resolve()),
        )
        self.assertEqual(
            command[command.index("--scaffolding_bam") + 1],
            str(bam.resolve()),
        )
        self.assertNotIn("--reassembly_dir", command)
        self.assertNotIn("--alignment_dir", command)

    def test_annotation_forwards_the_explicit_mag_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "annotation",
                    "--samplesheet",
                    str(root / "samples.csv"),
                    "--outdir",
                    str(root / "annotation-run"),
                    "--mag-dir",
                    str(root / "mags"),
                ]
            )
            command = manager.build_nextflow_command(args)
        self.assertEqual(
            command[command.index("--annotation_mag_dir") + 1],
            str((root / "mags").resolve()),
        )
        self.assertNotIn("--reassembly_dir", command)

    def test_numbered_entry_names_are_rejected(self) -> None:
        with self.assertRaises(argparse.ArgumentTypeError):
            manager.parse_entry_module("module3")

    def test_preprocessing_selection_is_configuration_only(self) -> None:
        help_text = manager.build_parser()._subparsers._group_actions[0].choices[
            "run"
        ].format_help()
        self.assertNotIn("--preprocessing-libraries", help_text)
        self.assertIn("-t THREADS", help_text)
        self.assertIn("-m SIZE", help_text)
        self.assertIn("--memory", help_text)

    def test_resource_override_aliases_are_validated(self) -> None:
        parsed = manager.build_parser().parse_args(
            [
                "run",
                "--samplesheet",
                "samples.csv",
                "--outdir",
                "results",
                "-t",
                "12",
                "-m",
                "96g",
            ]
        )
        self.assertEqual(parsed.threads, 12)
        self.assertEqual(parsed.memory, "96 GB")
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                manager.build_parser().parse_args(
                    [
                        "run",
                        "--samplesheet",
                        "samples.csv",
                        "--outdir",
                        "results",
                        "--threads",
                        "0",
                    ]
                )
            with self.assertRaises(SystemExit):
                manager.build_parser().parse_args(
                    [
                        "run",
                        "--samplesheet",
                        "samples.csv",
                        "--outdir",
                        "results",
                        "--memory",
                        "96",
                    ]
                )

    def test_selected_stage_requires_only_its_database(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample\nexample\n")
            gtdbtk = root / "gtdbtk"
            gtdbtk.mkdir()
            bins = root / "bins"
            bins.mkdir()
            (bins / "bin.1.fa").write_text(">contig\nACGT\n")
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "annotation",
                    "--samplesheet",
                    str(samplesheet),
                    "--outdir",
                    str(root / "results"),
                    "--gtdbtk-db",
                    str(gtdbtk),
                    "--mag-dir",
                    str(bins),
                ]
            )
            manager.validate_run_inputs(args)

    def test_annotation_requires_valid_sample_specific_mag_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample\ngut_01\ngut_02\n")
            gtdbtk = root / "gtdbtk"
            gtdbtk.mkdir()
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "annotation",
                    "--samplesheet",
                    str(samplesheet),
                    "--outdir",
                    str(root / "results"),
                    "--gtdbtk-db",
                    str(gtdbtk),
                ]
            )
            with self.assertRaisesRegex(manager.MetahictError, "--mag-dir"):
                manager.validate_run_inputs(args)

            shared = root / "bins"
            shared.mkdir()
            (shared / "bin.fa").write_text(">contig\nACGT\n")
            args.mag_dir = str(shared)
            with self.assertRaisesRegex(manager.MetahictError, r"requires \{sample\}"):
                manager.validate_run_inputs(args)

            template = root / "{sample}" / "reassembled_bins"
            args.mag_dir = str(template)
            with self.assertRaisesRegex(manager.MetahictError, "not found"):
                manager.validate_run_inputs(args)

            for sample in ("gut_01", "gut_02"):
                directory = root / sample / "reassembled_bins"
                directory.mkdir(parents=True)
            with self.assertRaisesRegex(manager.MetahictError, "contains no"):
                manager.validate_run_inputs(args)

            for sample in ("gut_01", "gut_02"):
                (root / sample / "reassembled_bins" / f"{sample}.fa").write_text(
                    ">contig\nACGT\n"
                )
            manager.validate_run_inputs(args)

    def test_mge_requires_explicit_hosts_and_hic_or_contacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample,enzyme\ngut_01,Sau3AI\n")
            assembly = root / "assembly.fa"
            assembly.write_text(">contig\nACGT\n")
            hosts = root / "hosts"
            hosts.mkdir()
            (hosts / "host_mag.fa").write_text(">contig\nACGT\n")
            genomad = root / "genomad"
            genomad.mkdir()
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "mge",
                    "--samplesheet",
                    str(samplesheet),
                    "--fasta",
                    str(assembly),
                    "--outdir",
                    str(root / "results"),
                    "--genomad-db",
                    str(genomad),
                ]
            )
            with self.assertRaisesRegex(manager.MetahictError, "--host-dir"):
                manager.validate_run_inputs(args)

            args.host_dir = str(hosts)
            with self.assertRaisesRegex(manager.MetahictError, "--preprocessing-dir"):
                manager.validate_run_inputs(args)

            preprocessing = root / "preprocessing"
            (preprocessing / "hic").mkdir(parents=True)
            args.preprocessing_dir = str(preprocessing)
            manager.validate_run_inputs(args)
            command = manager.build_nextflow_command(args)
        self.assertEqual(
            command[command.index("--mge_host_dir") + 1],
            str(hosts.resolve()),
        )
        self.assertEqual(
            command[command.index("--mge_fasta") + 1],
            str(assembly.resolve()),
        )

    def test_preprocessing_directory_validates_required_library_children(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample,enzyme\ngut_01,Sau3AI\n")
            preprocessing = root / "1_preprocessing"
            (preprocessing / "sg").mkdir(parents=True)
            args = argparse.Namespace(preprocessing_dir=str(preprocessing))

            with self.assertRaisesRegex(manager.MetahictError, r"1_preprocessing/hic"):
                manager.validate_preprocessing_input(args, samplesheet, "reassembly")

            (preprocessing / "hic").mkdir()
            manager.validate_preprocessing_input(args, samplesheet, "reassembly")

    def test_scaffolding_requires_an_existing_explicit_bin(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample,enzyme\ngut_01,Sau3AI\n")
            checkm2 = root / "checkm2.dmnd"
            checkm2.touch()
            preprocessing = root / "preprocessing"
            (preprocessing / "hic").mkdir(parents=True)
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "scaffolding",
                    "--samplesheet",
                    str(samplesheet),
                    "--outdir",
                    str(root / "results"),
                    "--checkm2-db",
                    str(checkm2),
                    "--preprocessing-dir",
                    str(preprocessing),
                ]
            )
            with self.assertRaisesRegex(manager.MetahictError, "--scaffolding-bin"):
                manager.validate_run_inputs(args)

            missing = root / "missing.fa"
            args.scaffolding_bin = str(missing)
            with self.assertRaisesRegex(manager.MetahictError, "not found"):
                manager.validate_run_inputs(args)

            missing.write_text(">contig\nACGT\n")
            manager.validate_run_inputs(args)

            args.scaffolding_bam = str(root / "missing.bam")
            with self.assertRaisesRegex(manager.MetahictError, "Scaffolding BAM not found"):
                manager.validate_run_inputs(args)
            args.scaffolding_bam = None

            args.alignment_dir = str(root / "old-alignment")
            with self.assertRaisesRegex(manager.MetahictError, "does not use --alignment-dir"):
                manager.validate_run_inputs(args)

    def test_scaffolding_skip_checkm2_removes_database_requirement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample,enzyme\ngut_01,Sau3AI\n")
            bin_fasta = root / "bin.fa"
            bin_fasta.write_text(">contig\nACGT\n")
            preprocessing = root / "preprocessing"
            (preprocessing / "hic").mkdir(parents=True)
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "scaffolding",
                    "--samplesheet",
                    str(samplesheet),
                    "--outdir",
                    str(root / "results"),
                    "--scaffolding-bin",
                    str(bin_fasta),
                    "--preprocessing-dir",
                    str(preprocessing),
                    "--scaffolding-skip-checkm2",
                ]
            )
            manager.validate_run_inputs(args)
            command = manager.build_nextflow_command(args)

        self.assertIn("--scaffolding_skip_checkm2", command)
        self.assertNotIn("--checkm2_db", command)

    def test_scientific_entries_reject_a_missing_enzyme(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample,enzyme\ngut_01,\n")
            bin_fasta = root / "bin.fa"
            bin_fasta.write_text(">contig\nACGT\n")
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "scaffolding",
                    "--samplesheet",
                    str(samplesheet),
                    "--outdir",
                    str(root / "results"),
                    "--scaffolding-bin",
                    str(bin_fasta),
                    "--scaffolding-skip-checkm2",
                ]
            )
            with self.assertRaisesRegex(manager.MetahictError, "Missing restriction enzyme"):
                manager.validate_run_inputs(args)

    def test_obsolete_workflow_configuration_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text("sample,sg1,sg2,hic1,hic2\ngut_01,a,b,c,d\n")
            config = root / "metahict_configuration.yaml"
            config.write_text(
                'workflow:\n  restriction_enzyme: "Sau3AI,MluCI"\nmodules: {}\n'
            )
            bin_fasta = root / "bin.fa"
            bin_fasta.write_text(">contig\nACGT\n")
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "scaffolding",
                    "--samplesheet",
                    str(samplesheet),
                    "--config",
                    str(config),
                    "--outdir",
                    str(root / "results"),
                    "--scaffolding-bin",
                    str(bin_fasta),
                    "--scaffolding-skip-checkm2",
                ]
            )
            with self.assertRaisesRegex(manager.MetahictError, "Obsolete configuration schema"):
                manager.validate_run_inputs(args)

            config.write_text(
                "resources:\n"
                '  mge: {threads: 16, memory: "32 GB"}\n'
                '  mge_alignment: {threads: 16, memory: "32 GB"}\n'
                '  mge_contact: {threads: 1, memory: "32 GB"}\n'
                "modules: {}\n"
            )
            with self.assertRaisesRegex(manager.MetahictError, "resources.mge_alignment"):
                manager.validate_configuration_schema(config)

            config.write_text(
                "resources:\n"
                '  mge: {threads: 16, memory: "32 GB"}\n'
                "modules:\n"
                "  mge_alignment:\n"
                "    sort_memory: 1G\n"
                "  mge_contact:\n"
                "    method: normcc\n"
            )
            with self.assertRaisesRegex(manager.MetahictError, "modules.mge_alignment"):
                manager.validate_configuration_schema(config)

            config.write_text(
                "resources:\n"
                "  contact:\n"
                "    threads: 1\n"
                '    memory: "32 GB"\n'
                "modules:\n"
                "  contact:\n"
                "    normcc:\n"
                "      epsilonn: 1\n"
            )
            with self.assertRaisesRegex(
                manager.MetahictError, "unknown key: modules.contact.normcc.epsilonn"
            ):
                manager.validate_configuration_schema(config)

    def test_container_profiles_are_not_public_run_options(self) -> None:
        with io.StringIO() as errors, contextlib.redirect_stderr(errors):
            with self.assertRaises(SystemExit):
                manager.build_parser().parse_args(
                    [
                        "run",
                        "--samplesheet",
                        "samples.csv",
                        "--outdir",
                        "results",
                        "--profile",
                        "docker",
                    ]
                )


class RunLoggingTest(unittest.TestCase):
    class FakeProcess:
        def __init__(self, status: int, output: str = "Nextflow task output\n") -> None:
            self.status = status
            self.stdout = io.StringIO(output)

        def wait(self) -> int:
            return self.status

    def run_args(self, root: Path) -> argparse.Namespace:
        samplesheet = root / "samples.csv"
        samplesheet.write_text("sample,sg1,sg2,hic1,hic2\ngut_01,a,b,c,d\n")
        return manager.build_parser().parse_args(
            [
                "run",
                "--entry-module",
                "preprocessing",
                "--samplesheet",
                str(samplesheet),
                "--outdir",
                str(root / "results"),
                "--report-dir",
                str(root / "reports"),
            ]
        )

    def test_run_preflight_is_concise_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = self.run_args(Path(temporary))
            with mock.patch.object(manager, "verify_runtime") as verify:
                manager.validate_run_inputs(args)
        verify.assert_called_once_with(verbose=False)

    def test_verbose_preflight_is_forwarded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = self.run_args(Path(temporary))
            args.verbose_preflight = True
            with mock.patch.object(manager, "verify_runtime") as verify:
                manager.validate_run_inputs(args)
        verify.assert_called_once_with(verbose=True)

    def test_external_preflight_output_is_tee_logged(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            log_path = Path(temporary) / "run.log"
            terminal_stdout = io.StringIO()
            terminal_stderr = io.StringIO()
            with log_path.open("w", buffering=1) as log_handle, contextlib.redirect_stdout(
                manager.TeeTextIO(terminal_stdout, log_handle)
            ), contextlib.redirect_stderr(manager.TeeTextIO(terminal_stderr, log_handle)):
                manager.run_command(
                    [
                        sys.executable,
                        "-c",
                        "import sys; print('preflight stdout'); print('preflight stderr', file=sys.stderr)",
                    ]
                )
            logged = log_path.read_text()
            self.assertIn("preflight stdout", logged)
            self.assertIn("preflight stderr", logged)

    def test_run_directories_do_not_collide_or_delete_legacy_log(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            reports = Path(temporary) / "reports"
            reports.mkdir()
            (reports / "run.log").write_text("legacy log\n")
            fixed = manager.datetime(2026, 8, 27, 14, 30, 15, tzinfo=manager.timezone.utc)
            first_id, first = manager.create_run_directory(reports, fixed)
            manager.prepare_report_links(reports, first_id, first)
            second_id, second = manager.create_run_directory(reports, fixed)
            manager.prepare_report_links(reports, second_id, second)

            self.assertEqual((first_id, second_id), ("20260827T143015Z", "20260827T143015Z-001"))
            self.assertTrue(first.is_dir())
            self.assertTrue(second.is_dir())
            self.assertEqual((reports / "latest").resolve(), second.resolve())
            legacy_logs = list((reports / "runs").glob("legacy-before-*/run.log"))
            self.assertEqual(len(legacy_logs), 1)
            self.assertEqual(legacy_logs[0].read_text(), "legacy log\n")

    def test_completed_run_is_immutable_and_captures_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.run_args(root)
            stdout = io.StringIO()
            with mock.patch.object(
                manager, "validate_run_inputs", side_effect=lambda _: print("[PASS] preflight")
            ), mock.patch.object(manager, "git_revision", return_value="abc123"), mock.patch.object(
                manager.platform, "platform", return_value="test-platform"
            ), mock.patch.object(
                manager.subprocess, "Popen", return_value=self.FakeProcess(0)
            ), contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
                manager.command_run(args)

            reports = root / "reports"
            run_directories = list((reports / "runs").glob("20*"))
            self.assertEqual(len(run_directories), 1)
            run_directory = run_directories[0]
            run_log = (run_directory / "run.log").read_text()
            self.assertIn("[PASS] preflight", run_log)
            self.assertIn("Nextflow task output", run_log)
            self.assertIn("[PASS] METAHICT workflow completed", run_log)
            self.assertTrue((reports / "run.log").is_symlink())
            self.assertEqual((reports / "latest").resolve(), run_directory.resolve())
            provenance = json.loads((run_directory / "provenance.json").read_text())
            self.assertEqual((provenance["status"], provenance["exit_code"]), ("completed", 0))
            self.assertEqual(provenance["git_revision"], "abc123")
            self.assertTrue((run_directory / "software_versions.tsv").is_file())
            self.assertIn(
                "Binning-refiner\t1.4.3\twheel",
                (run_directory / "software_versions.tsv").read_text(),
            )
            self.assertTrue((run_directory / "process_logs" / "index.tsv").is_file())

    def test_failed_run_retains_log_and_failure_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.run_args(root)
            with mock.patch.object(manager, "validate_run_inputs"), mock.patch.object(
                manager, "git_revision", return_value=None
            ), mock.patch.object(
                manager.platform, "platform", return_value="test-platform"
            ), mock.patch.object(
                manager.subprocess,
                "Popen",
                return_value=self.FakeProcess(7, "task failed\n"),
            ), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(manager.MetahictError):
                    manager.command_run(args)

            run_directory = next((root / "reports" / "runs").glob("20*"))
            provenance = json.loads((run_directory / "provenance.json").read_text())
            self.assertEqual((provenance["status"], provenance["exit_code"]), ("failed", 7))
            self.assertIn("task failed", (run_directory / "run.log").read_text())
            self.assertIn("Exit code: 7", (run_directory / "failure_summary.txt").read_text())

    def test_missing_samplesheet_is_recorded_as_preflight_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = self.run_args(root)
            Path(args.samplesheet).unlink()
            with mock.patch.object(manager, "git_revision", return_value=None), mock.patch.object(
                manager.platform, "platform", return_value="test-platform"
            ), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(manager.MetahictError):
                    manager.command_run(args)

            run_directory = next((root / "reports" / "runs").glob("20*"))
            self.assertTrue((run_directory / "samplesheet.sha256").read_text().startswith("MISSING  "))
            self.assertIn("Samplesheet not found", (run_directory / "run.log").read_text())
            provenance = json.loads((run_directory / "provenance.json").read_text())
            self.assertEqual(provenance["status"], "failed")
            self.assertIsNone(provenance["samplesheet_sha256"])

    def test_process_logs_are_archived_from_trace_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run_directory = root / "run"
            work = root / "work"
            run_directory.mkdir()
            (run_directory / "trace.txt").write_text(
                "task_id\thash\tname\tstatus\texit\n"
                "1\t9d/68d62c\tPREPROCESSING_SG (gut_01:sg)\tCOMPLETED\t0\n"
            )
            task = work / "9d" / "68d62cabcdef"
            task.mkdir(parents=True)
            (task / ".command.out").write_text("tool stdout\n")
            (task / ".command.err").write_text("tool stderr\n")
            (task / ".command.sh").write_text("python preprocessing.py\n")
            (task / ".exitcode").write_text("0\n")

            archived = manager.archive_process_logs(run_directory, work)

            self.assertEqual(len(archived), 1)
            log_directory = Path(archived[0]["log_directory"])
            self.assertEqual((log_directory / "command.out").read_text(), "tool stdout\n")
            self.assertEqual((log_directory / "command.err").read_text(), "tool stderr\n")
            self.assertIn("PREPROCESSING_SG", (run_directory / "process_logs" / "index.tsv").read_text())


class SamplesheetTest(unittest.TestCase):
    def test_samplesheet_command_writes_absolute_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = []
            for name in ("sg_R1.fastq.gz", "sg_R2.fastq.gz", "hic_R1.fastq.gz", "hic_R2.fastq.gz"):
                path = root / name
                path.touch()
                inputs.append(path)
            output = root / "samples.csv"
            args = argparse.Namespace(
                sample="gut_01",
                sg_r1=str(inputs[0]),
                sg_r2=str(inputs[1]),
                hic_r1=str(inputs[2]),
                hic_r2=str(inputs[3]),
                enzyme="Sau3AI,MluCI",
                output=str(output),
                allow_missing=False,
                append=False,
            )
            manager.command_samplesheet(args)
            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle))
        self.assertEqual(row["sample"], "gut_01")
        self.assertEqual(
            list(row),
            ["sample", "sg1", "sg2", "hic1", "hic2", "enzyme", "long_read_type"],
        )
        self.assertEqual(row["enzyme"], "Sau3AI,MluCI")
        self.assertEqual(row["long_read_type"], "")
        self.assertTrue(Path(row["sg1"]).is_absolute())

    def test_samplesheet_command_writes_single_long_read_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            long_reads = root / "nanopore.fastq.gz"
            hic1 = root / "hic_R1.fastq.gz"
            hic2 = root / "hic_R2.fastq.gz"
            for path in (long_reads, hic1, hic2):
                path.touch()
            output = root / "samples.csv"
            manager.command_samplesheet(
                argparse.Namespace(
                    sample="long_sample",
                    sg_r1=str(long_reads),
                    sg_r2=None,
                    long_read_type="nano-hq",
                    hic_r1=str(hic1),
                    hic_r2=str(hic2),
                    enzyme="Sau3AI",
                    output=str(output),
                    allow_missing=False,
                    append=False,
                )
            )
            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle))
        self.assertEqual(row["sg2"], "")
        self.assertEqual(row["long_read_type"], "nano-hq")

    def test_long_read_assembly_does_not_require_preprocessing_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.csv"
            samplesheet.write_text(
                "sample,sg1,sg2,hic1,hic2,enzyme,long_read_type\n"
                "long_sample,/reads/long.fastq.gz,,/reads/hic1.fastq.gz,"
                "/reads/hic2.fastq.gz,Sau3AI,nano-hq\n"
            )
            args = manager.build_parser().parse_args(
                [
                    "run",
                    "--entry-module",
                    "assembly",
                    "--samplesheet",
                    str(samplesheet),
                    "--config",
                    str(manager.DEFAULT_CONFIG),
                    "--outdir",
                    str(root / "results"),
                ]
            )
            with mock.patch.object(manager, "verify_runtime"):
                manager.validate_run_inputs(args)
            command = manager.build_nextflow_command(args)
        self.assertNotIn("--preprocessing_dir", command)

    def test_config_command_uses_the_standard_filename_and_has_no_enzyme(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "metahict_configuration.yaml"
            manager.command_config(
                argparse.Namespace(output=str(output), force=False)
            )
            text = output.read_text()
        self.assertNotIn("restriction_enzyme", text)
        self.assertNotIn("workflow:", text)
        self.assertNotIn("{cpus:", text)
        self.assertNotIn("time:", text)
        self.assertIn("preprocessing: {threads: 8, memory: \"32 GB\"}", text)
        self.assertIn("contact: {threads: 1, memory: \"32 GB\"}", text)
        self.assertIn("scaffolding: {threads: 8, memory: \"32 GB\"}", text)
        self.assertIn("annotation: {threads: 8, memory: \"64 GB\"}", text)
        resources = text.split("resources:\n", 1)[1].split("\nmodules:\n", 1)[0]
        self.assertIn("mge: {threads: 16, memory: \"32 GB\"}", resources)
        self.assertNotIn("mge_alignment:", resources)
        self.assertNotIn("mge_contact:", resources)
        self.assertNotIn("\n  mge_alignment:", text)
        self.assertNotIn("\n  mge_contact:", text)
        self.assertIn("\n    alignment:\n", text)
        self.assertIn("\n    contact:\n", text)
        self.assertIn("      shotgun:\n        enabled: true\n        deduplicate: false", text)
        self.assertIn("      hic:\n        enabled: true\n        deduplicate: true", text)


class CompareCommandTest(unittest.TestCase):
    def test_compare_parser_has_release_manifest_default(self) -> None:
        args = manager.build_parser().parse_args(
            ["compare", "--baseline", "old", "--candidate", "new"]
        )
        self.assertIsNone(args.manifest)
        self.assertEqual(args.report, "validation/scientific-regression.tsv")


class StubInputTest(unittest.TestCase):
    def test_stub_fixture_is_self_contained_and_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "stub_inputs"
            result = subprocess.run(
                [
                    sys.executable,
                    str(manager.PROJECT_ROOT / "nextflow" / "ci" / "create_stub_inputs.py"),
                    "--output-dir",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            samplesheet = Path(result.stdout.strip())
            with samplesheet.open(newline="") as handle:
                rows = list(csv.DictReader(handle))

            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["sample"], "example_dataset")
            self.assertEqual(rows[0]["enzyme"], "Sau3AI,MluCI")
            for column in ("sg1", "sg2", "hic1", "hic2"):
                read_file = Path(rows[0][column])
                self.assertTrue(read_file.is_absolute())
                with gzip.open(read_file, "rt") as handle:
                    self.assertTrue(handle.readline().startswith("@stub_read_"))
            self.assertEqual(rows[0]["long_read_type"], "")
            self.assertEqual(rows[1]["sample"], "long_read_example")
            self.assertEqual(rows[1]["sg2"], "")
            self.assertEqual(rows[1]["long_read_type"], "nano-hq")
            with gzip.open(Path(rows[1]["sg1"]), "rt") as handle:
                self.assertTrue(handle.readline().startswith("@stub_long_read"))


class ExampleScaffoldingTest(unittest.TestCase):
    def test_workflow_scope_runs_the_dependency_free_stub_profile(self) -> None:
        args = argparse.Namespace(scope="workflow")
        with mock.patch.object(manager, "run_stub_test") as run_stub:
            manager.command_test(args)
        run_stub.assert_called_once_with("stub")

    def test_discovers_every_supported_bin_without_assuming_a_count(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            short_dir = (
                results
                / "short_sample"
                / "7_reassembly"
                / "reassembled_bins"
            )
            short_dir.mkdir(parents=True)
            for name in ("bin1.fa", "bin2.fasta", "bin3.fna", "notes.txt"):
                (short_dir / name).touch()
            bins = manager.example_scaffolding_bins(
                results, {"sample": "short_sample", "long_read_type": ""}
            )
        self.assertEqual(
            [path.name for path in bins], ["bin1.fa", "bin2.fasta", "bin3.fna"]
        )

    def test_dynamic_validation_accepts_completed_and_skipped_bins(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            completed_bin = results / "inputs" / "completed.fa"
            skipped_bin = results / "inputs" / "skipped.fa"
            completed_bin.parent.mkdir()
            completed_bin.touch()
            skipped_bin.touch()
            for bin_fasta, status in (
                (completed_bin, "completed"),
                (skipped_bin, "skipped"),
            ):
                output = results / "sample" / "8_scaffolding" / bin_fasta.stem
                output.mkdir(parents=True)
                (output / "scaffolding_status.tsv").write_text(
                    "bin\tstatus\treason\n"
                    f"{bin_fasta.name}\t{status}\t\n"
                )
            manager.validate_example_scaffolding(
                results,
                [("sample", completed_bin), ("sample", skipped_bin)],
            )

    def test_example_scope_has_simple_defaults(self) -> None:
        args = manager.build_parser().parse_args(["test", "example"])
        self.assertEqual(args.scope, "example")
        self.assertEqual(args.outdir, "results")
        self.assertTrue(str(args.config).endswith("example_dataset_configuration.yaml"))


class TopLevelHelpTest(unittest.TestCase):
    def test_top_level_help_is_a_first_run_guide(self) -> None:
        help_text = manager.build_parser().format_help()
        self.assertIn("usage: ./metahict", help_text)
        self.assertIn("commands:", help_text)
        self.assertIn("GETTING STARTED", help_text)
        self.assertIn("./metahict install", help_text)
        self.assertIn("./metahict database all", help_text)
        self.assertIn("./metahict doctor --runtime --databases", help_text)
        self.assertIn("./metahict run --help", help_text)
        self.assertIn("./metahict run --entry-module MODULE --help", help_text)


class GeneralRunHelpTest(unittest.TestCase):
    def test_run_help_explains_configuration_resources_and_module_help(self) -> None:
        run_parser = manager.build_parser()._subparsers._group_actions[0].choices[
            "run"
        ]
        help_text = run_parser.format_help()
        self.assertIn("required and workflow selection", help_text)
        self.assertIn("run-wide resource overrides", help_text)
        self.assertNotIn("--profile", help_text)
        self.assertNotIn("--container", help_text)
        self.assertIn("reference-database overrides", help_text)
        self.assertIn("reusable upstream results for selected-module runs", help_text)
        self.assertNotIn("--input-root", help_text)
        self.assertIn("Read-file paths and restriction enzymes come from --samplesheet", help_text)
        self.assertIn("Algorithm settings and normal per-module resources come from --config", help_text)
        self.assertIn("values override the YAML resources for one run", help_text)
        self.assertIn("./metahict run --entry-module binning --help", help_text)
        self.assertIn("--preprocessing-dir DIR", help_text)
        self.assertNotIn("--sg-preprocessing-dir", help_text)
        self.assertNotIn("--hic-preprocessing-dir", help_text)
        self.assertIn("--assembly-dir DIR", help_text)
        self.assertIn("--mag-dir DIR", help_text)
        self.assertIn("--fasta FASTA", help_text)
        self.assertIn("--host-dir DIR", help_text)


class ModuleHelpTest(unittest.TestCase):
    def test_every_entry_has_detailed_help_without_required_run_inputs(self) -> None:
        expected_keys = {
            "preprocessing": "modules.preprocessing.libraries.hic.deduplicate",
            "assembly": "modules.assembly.assembler",
            "alignment": "modules.alignment.bwa.options",
            "coverage": "modules.coverage.mapping.min_percent_identity",
            "contact": "modules.contact.method",
            "binning": "modules.binning.refinement.min_completeness",
            "reassembly": "modules.reassembly.read_selection.em_cutoff_quantile",
            "scaffolding": "modules.scaffolding.yahs.min_mapping_quality",
            "annotation": "modules.annotation.pplacer_threads",
            "mge": "modules.mge.genomad.splits",
        }
        for module, expected_key in expected_keys.items():
            with self.subTest(module=module), io.StringIO() as output:
                with contextlib.redirect_stdout(output):
                    status = manager.main(
                        ["run", "--entry-module", module, "--help"]
                    )
                help_text = output.getvalue()
                self.assertEqual(status, 0)
                self.assertIn("MODULE-SPECIFIC HELP", help_text)
                self.assertIn("RESOURCE BEHAVIOR", help_text)
                self.assertIn(f"resources.{module}.threads=", help_text)
                self.assertIn(expected_key, help_text)
                for key in manager.configured_module_keys(module):
                    self.assertIn(f"modules.{module}.{key}", help_text)
                self.assertIn("RECOMMENDED COMMAND", help_text)

    def test_short_help_alias_selects_the_same_module_help(self) -> None:
        with io.StringIO() as output, contextlib.redirect_stdout(output):
            status = manager.main(
                ["run", "--entry-module=annotation", "-h"]
            )
            help_text = output.getvalue()
        self.assertEqual(status, 0)
        self.assertIn("resources.annotation.memory=64 GB", help_text)
        self.assertIn("--mag-dir", help_text)


if __name__ == "__main__":
    unittest.main()
