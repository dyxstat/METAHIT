# Module documentation

METAHICT can be run as a full workflow or as selected module-level entry workflows. The module pages below summarize each module’s purpose, inputs, outputs, parameters, and example commands.

Each module page shows two command styles:

- The Nextflow command is the recommended interface for routine runs. It uses the sample sheet, selected execution profile, standardized output layout, and Nextflow provenance reports.
- The manual wrapper command calls `metahict.py` directly. 

| Module | Page |
| --- | --- |
| Module 1 | [Preprocessing](module1_preprocessing.md) |
| Module 2 | [Assembly](module2_assembly.md) |
| Module 3 | [Alignment](module3_alignment.md) |
| Module 4 | [Coverage](module4_coverage.md) |
| Module 5 | [Contact generation and normalization](module5_contact.md) |
| Module 6 | [Binning and consolidation](module6_binning.md) |
| Module 7 | [Reassembly](module7_reassembly.md) |
| Module 8 | [Scaffolding](module8_scaffolding.md) |
| Module 9 | [Annotation](module9_annotation.md) |
| Module 10 | [MGE analysis](module10_mge.md) |

For the complete defaults file, see:

```text
nextflow/assets/default_params.yaml
```

For full wrapper help, run:

```bash
python metahict.py <module> --help
```

For contact normalization, run:

```bash
python metahict.py contact normcc --help
```
