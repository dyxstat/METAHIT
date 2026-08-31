# Contact generation and normalization

The contact module converts the filtered Hi-C BAM into a contig-by-contig
matrix, removes weak or potentially spurious signals, and applies the selected
normalization. NormCC is the default method.

Raw contact counts are affected by contig length, restriction sites, coverage,
mapping ambiguity, and experimental noise. Normalization makes contacts more
comparable. Genome membership still requires evidence from binning and quality
assessment.

## Before you run

The BAM must have been aligned to the exact input assembly. The samplesheet
must contain the restriction enzyme or comma-separated enzyme list used for
this Hi-C library. HiCzin and MetaTOR coverage files must come from the same
assembly and sample.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--assembly-dir` | `./metahict run` | Assembly output directory containing the contig FASTA |
| `--alignment-dir` | `./metahict run` | Hi-C alignment output directory |
| `enzyme` | Samplesheet | Required restriction enzyme list |

## Optional parameters

Set these stage options under `modules.contact` in `metahict_configuration.yaml`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `method` | `normcc` | Contact-normalization method selected by the Nextflow defaults |
| `raw_contacts.min_contact_signal` | `1` | Minimum contact signal during raw-matrix filtering |
| `raw_contacts.min_contig_length` | `1000` | Minimum contig length |
| `raw_contacts.min_mapping_quality` | `30` | Minimum alignment MAPQ |
| `raw_contacts.min_aligned_length` | `30` | Minimum aligned match length |
| `denoising.spurious_contact_percent` | `5` | Lowest normalized-contact percentage removed as potentially spurious |
| `normcc.epsilon` | `1` | Numerical stabilizer used by NormCC |
| `hiczin.coverage_file` | - | Absolute path to the matching coverage table required when `method: hiczin` |
| `hiczin.epsilon` | `1` | Numerical stabilizer used by HiCzin |
| `bin3c.epsilon` | `1` | Numerical stabilizer used by bin3C normalization |
| `bin3c.max_iterations` | `1000` | Maximum bin3C balancing iterations |
| `bin3c.convergence_tolerance` | `1e-6` | bin3C balancing convergence tolerance |
| `metator.coverage_file` | - | Absolute path to the matching coverage table required when `method: metator` |
| `metator.epsilon` | `1` | Numerical stabilizer used by MetaTOR normalization |

Valid methods are `raw`, `normcc`, `hiczin`, `bin3c`, and `metator`. Only the
subsection for the selected method is passed. `hiczin` and `metator` require a
coverage table whose contig identifiers match the input FASTA.

`raw_contacts` settings are applied while constructing the observed matrix.
`denoising` is applied after normalization. Only the map named by `method` is
active; changing `hiczin.epsilon` has no effect when `method: normcc`.

## Outputs

| Output | Description |
| --- | --- |
| `5_contact/Raw_contact_matrix.npz` | Raw contig-by-contig contact counts |
| `5_contact/denoised_contact_matrix_<METHOD>.npz` | Matrix produced by the selected normalization method |

## How to check the result

Confirm that the raw and selected-method matrices are non-empty. Compare the
number of retained contigs and contacts before and after filtering. A nearly
empty matrix can result from a wrong enzyme, mismatched BAM/assembly, low Hi-C
mapping, short/fragmented contigs, or overly strict thresholds.

## Recommended command

```bash
./metahict run \
  --entry-module contact \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --assembly-dir results/sample_01/2_assembly \
  --alignment-dir results/sample_01/3_alignment \
  --outdir results/contact
```

[Previous: Coverage](coverage.md) · [Next: Binning](binning.md) · [All modules](README.md)
