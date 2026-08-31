# Third-party software provenance

This manifest records the direct external programs invoked by the current
METAHICT workflow. METAHICT does not store compiled third-party executables;
it installs programs from their original distributors, Bioconda/Conda-Forge
package artifacts, or the attributed GitHub forks below. Source-level code
adaptations maintained in METAHICT are identified separately in the notes.

The canonical, complete record of **every resolved package**, including Python
libraries and operating-system dependencies, is in the explicit files under
`installation/locks/`. Each line is an exact channel/package/version/build
artifact. `installation/locks/SHA256SUMS`
checks that those lock files have not changed.  This table records the
scientific tools that METAHICT directly calls and therefore need a source,
exact version or commit, licence, citation, and installation route.

| Software | Upstream source | Exact METAHICT release input | Licence | Citation | Installation route |
| --- | --- | --- | --- | --- | --- |
| Nextflow | https://github.com/nextflow-io/nextflow | 26.04.6 | Apache-2.0 | Di Tommaso P *et al.* 2017. https://doi.org/10.1038/nbt.3820 | Version-pinned launcher in `nextflow/bin/nextflow`; Java runtime in the locked `metahict_nextflow_env` |
| FastQC | https://github.com/s-andrews/FastQC | `fastqc` 0.12.1, Bioconda build `hdfd78af_0` | GPL-3.0 | Andrews S. *FastQC*; https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ | Locked Conda package (`metahict_env.explicit.txt`) |
| BBTools / BBDuk | https://sourceforge.net/projects/bbmap/ | `bbmap` 39.10, Bioconda build `h92535d8_0` | BSD 3-Clause | Bushnell B. *BBTools*; https://bbmap.org/ | Locked Conda package (`metahict_env.explicit.txt`) |
| Seqtk | https://github.com/lh3/seqtk | `seqtk` 1.4, Bioconda build `he4a0461_2` | MIT | Li H. *Seqtk* software repository | Locked Conda package (`metahict_env.explicit.txt`) |
| BWA | https://github.com/lh3/bwa | `bwa` 0.7.18, Bioconda build `he4a0461_1` | GPL-3.0 | Li H, Durbin R. 2009. https://doi.org/10.1093/bioinformatics/btp698 | Locked Conda package (`metahict_env.explicit.txt`) |
| SAMtools / HTSlib | https://github.com/samtools/samtools | `samtools` 1.21, Bioconda build `h50ea8bc_0` | MIT | Danecek P *et al.* 2021. https://doi.org/10.1093/gigascience/giab008 | Locked Conda package (`metahict_env.explicit.txt`) |
| MetaBAT2 (`jgi_summarize_bam_contig_depths`) | https://bitbucket.org/berkeleylab/metabat/ | `metabat2` 2.18, Bioconda build `h38e344b_2` | BSD-3-Clause | Kang DD *et al.* 2019. https://doi.org/10.7717/peerj.7359 | Locked Conda package (`metabat2.explicit.txt`); no executable is stored in `modules/` |
| MEGAHIT | https://github.com/voutcn/megahit | `megahit` 1.2.9, Bioconda build `h43eeafb_5` | GPL-3.0 | Li D *et al.* 2015. https://doi.org/10.1093/bioinformatics/btv033 | Locked Conda package (`metahict_env.explicit.txt`) |
| metaSPAdes | https://github.com/ablab/spades | `spades` 4.0.0, Bioconda build `h6dccd9a_3` | GPL-2.0-only | Bankevich A *et al.* 2012. https://doi.org/10.1089/cmb.2012.0021 | Locked Conda package (`metahict_env.explicit.txt`) |
| metaFlye | https://github.com/mikolmogorov/Flye | `flye` 2.9.5, Bioconda build `py39hdf45acc_1` | BSD-3-Clause | Kolmogorov M *et al.* 2019. https://doi.org/10.1038/s41587-019-0072-8 | Locked Conda package (`metahict_env.explicit.txt`) |
| QUAST | https://github.com/ablab/quast | `quast` 5.3.0, Bioconda build `py39pl5321h746d604_1` | Custom; see upstream distribution | Gurevich A *et al.* 2013. https://doi.org/10.1093/bioinformatics/btt086 | Locked Conda package (`metahict_env.explicit.txt`) |
| CheckM | https://github.com/Ecogenomics/CheckM | `checkm-genome` 1.2.4, Bioconda build `pyhdfd78af_0` | GPL-3.0-or-later | Parks DH *et al.* 2015. https://doi.org/10.1101/gr.186072.114 | Locked Conda package (`metahict_env.explicit.txt`) |
| CheckM2 | https://github.com/chklovski/CheckM2 | `checkm2` 1.1.0, Bioconda build `pyh7e72e81_1` | GPL-3.0-only | Chklovski A *et al.* 2023. https://doi.org/10.1038/s41592-023-01940-w | Locked Conda package (`checkm2.explicit.txt`) |
| FragGeneScan | https://sourceforge.net/projects/fraggenescan/ | `fraggenescan` 1.32, Bioconda build `h7b50bb2_1` | BSD | Rho M *et al.* 2010. https://doi.org/10.1093/nar/gkq747 | Locked Conda package (`metahict_env.explicit.txt`) |
| HMMER | http://hmmer.org/ | `hmmer` 3.4, Bioconda build `hb6cb901_4` | BSD | Eddy SR. 2011. https://doi.org/10.1371/journal.pcbi.1002195 | Locked Conda package (`metahict_env.explicit.txt`) |
| GTDB-Tk | https://github.com/Ecogenomics/GTDBTk | `gtdbtk` 2.4.0, Bioconda build `pyhdfd78af_2` | GPL-3.0 | Chaumeil PA *et al.* 2022. https://doi.org/10.1093/bioinformatics/btac672 | Locked Conda package (`gtdbtk-2.4.0.explicit.txt`) |
| geNomad | https://github.com/apcamargo/genomad | `genomad` 1.12.0, Bioconda build `pyhdfd78af_0` | BSD-4-Clause | Camargo AP *et al.* 2023. https://doi.org/10.1038/s41587-023-01953-y | Locked Conda package (`genomad.explicit.txt`) |
| ccfind | https://github.com/yosuken/ccfind | `ccfind` v1.4.7, Git commit `674366b49dd31cb909c2e52834e4ec8ede8919e7`; source-archive SHA-256 `abee22a03ab2e7c475d08361f90fa409a5626733649e4eb39b0cdbc7036ff463`; runtime dependencies in `ccfind_env` include `ruby` 3.4.8, `fasta3` 36.3.8i, BLAST 2.16.0, Prodigal 2.6.3, and GNU parallel 20260722 | MIT | Nishimura Y *et al.* 2017. https://doi.org/10.1128/mSphere.00359-16 | Verified source archive downloaded by `./metahict install`; runtime packages locked in `ccfind_env.explicit.txt` |
| DIAMOND | https://github.com/bbuchfink/diamond | 2.1.11 in `checkm2` | GPL-3.0-or-later | Buchfink B *et al.* 2021. https://doi.org/10.1038/s41592-021-01101-x | Locked Conda package (`checkm2.explicit.txt`) |
| Prodigal | https://github.com/hyattpd/Prodigal | 2.6.3 | GPL-3.0-or-later | Hyatt D *et al.* 2010. https://doi.org/10.1186/1471-2105-11-119 | Locked Conda package (`metahict_env.explicit.txt`); also present in the CheckM2, GTDB-Tk, and ccfind locks |
| YaHS | https://github.com/c-zhou/yahs | `yahs` 1.2.2, Bioconda build `he4a0461_0` | MIT | Zhou C *et al.* 2023. https://doi.org/10.1093/bioinformatics/btac808 | Locked Conda package (`metahict_env.explicit.txt`) |
| Binning_refiner | https://github.com/songweizhi/Binning_refiner | 1.4.3 PyPI wheel; SHA-256 `f5cb191492ff28298c9f2e9b6ed5299cebc2b33cbcc59f4e7ad3ebce7ab6c8a0` | AGPL-3.0-or-later | Song W, Thomas T. 2017. https://doi.org/10.1093/bioinformatics/btx086 | Installed from the upstream wheel pinned in `installation/pip-requirements.txt`; `run_binning_refiner.py` is a METAHICT output-layout adapter, not a copy of the algorithm |
| bin3C Python 3 port | https://github.com/dyxstat/bin3C branch `bin3C-Python3` (fork of https://github.com/cerebis/bin3C) | Git commit `181d80fb3f722d165a66cd3f782c578996e27675` | AGPL-3.0 | DeMaere MZ, Darling AE. 2019. https://doi.org/10.1186/s13059-019-1643-1 | Installed by `pip` from `installation/pip-requirements.txt` |
| MetaCC scripts for METAHICT | https://github.com/dyxstat/MetaCC branch `MetaCC-METAHICT` | Git commit `026a18f79f6a813625a8215855ef52e1b1a6ce7f` | See upstream repository | MetaCC software repository | Installed by `pip` from `installation/pip-requirements.txt` |
| ImputeCC scripts for METAHICT | https://github.com/dyxstat/ImputeCC branch `ImputeCC-METAHICT` | Git commit `9ad1119978d4e269a48106740ecd4b2882906cc1` | See upstream repository | ImputeCC software repository | Installed by `pip` from `installation/pip-requirements.txt` |

## Notes on METAHICT-owned code and data

Large reference databases are not distributed in the repository. Users obtain
them with `./metahict database` and provide their paths at run time. The
explicitly selected database releases include GTDB release 220. The geNomad
and CheckM2 download tools
record the database version in the downloaded database directory.

The bin3C fork retains upstream Git history and its AGPL-3.0 licence.  The
METAHICT environment installs bin3C, MetaCC, and ImputeCC from the immutable
commits above through `installation/pip-requirements.txt`; no bin3C, MetaCC,
or ImputeCC source tree is included in the active METAHICT module directories.

The bin-selection and argument-handling portions of
`modules/6_binning/bin_refinement.py` and
`modules/7_reassembly/reassemble_bins.py`, together with their small
helper scripts, were adapted from metaWRAP (GPL-3.0; Uritskiy GV *et al.*,
2018, https://doi.org/10.1186/s40168-018-0541-1) and are maintained as part of
METAHICT under its GPL license. They do not include the MetaBAT2 or
Binning_refiner programs. Unused duplicate metaWRAP helpers were removed.
