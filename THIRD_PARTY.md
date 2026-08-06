# Third-party software provenance

This manifest records the direct external programs invoked by the current
METAHICT workflow.  It is a record of provenance, not a copy of third-party
source code.  METAHICT installs the programs from their original distributors,
Bioconda/Conda-Forge package artifacts, or the attributed bin3C fork below.

The canonical, complete record of **every resolved package**, including Python
libraries and operating-system dependencies, is in
`installation/locks/linux-64/*.explicit.txt`.  Each line is an exact
channel/package/version/build artifact.  `installation/locks/SHA256SUMS`
checks that those lock files have not changed.  This table records the
scientific tools that METAHICT directly calls and therefore need a source,
licence, and citation.

| Software | Upstream source | Exact METAHICT release input | Licence | Citation |
| --- | --- | --- | --- | --- |
| BBTools / BBDuk | https://sourceforge.net/projects/bbmap/ | BBMap archive requested as `BBMap_39.10.tar.gz`; SHA-256 `ab5dfc0bbaa5be338596aec3558c7a7c891e8d8b186e9bd671552466215b9b15` | BSD 3-Clause | Bushnell B. *BBTools*; https://bbmap.org/ |
| BWA | https://github.com/lh3/bwa | `bwa` 0.7.18, Bioconda build `he4a0461_1` | GPL-3.0 | Li H, Durbin R. 2009. https://doi.org/10.1093/bioinformatics/btp698 |
| SAMtools / HTSlib | https://github.com/samtools/samtools | `samtools` 1.21, Bioconda build `h50ea8bc_0` | MIT | Danecek P *et al.* 2021. https://doi.org/10.1093/gigascience/giab008 |
| MEGAHIT | https://github.com/voutcn/megahit | `megahit` 1.2.9, Bioconda build `h43eeafb_5` | GPL-3.0 | Li D *et al.* 2015. https://doi.org/10.1093/bioinformatics/btv033 |
| metaSPAdes | https://github.com/ablab/spades | `spades` 4.0.0, Bioconda build `h6dccd9a_3` | GPL-2.0-only | Bankevich A *et al.* 2012. https://doi.org/10.1089/cmb.2012.0021 |
| metaFlye | https://github.com/mikolmogorov/Flye | `flye` 2.9.5, Bioconda build `py39hdf45acc_1` | BSD-3-Clause | Kolmogorov M *et al.* 2019. https://doi.org/10.1038/s41587-019-0072-8 |
| QUAST | https://github.com/ablab/quast | `quast` 5.3.0, Bioconda build `py39pl5321h746d604_1` | Custom; see upstream distribution | Gurevich A *et al.* 2013. https://doi.org/10.1093/bioinformatics/btt086 |
| CheckM | https://github.com/Ecogenomics/CheckM | `checkm-genome` 1.2.4, Bioconda build `pyhdfd78af_0` | GPL-3.0-or-later | Parks DH *et al.* 2015. https://doi.org/10.1101/gr.186072.114 |
| CheckM2 | https://github.com/chklovski/CheckM2 | `checkm2` 1.1.0, Bioconda build `pyh7e72e81_1` | GPL-3.0-only | Chklovski A *et al.* 2023. https://doi.org/10.1038/s41592-023-01940-w |
| GTDB-Tk | https://github.com/Ecogenomics/GTDBTk | `gtdbtk` 2.4.0, Bioconda build `pyhdfd78af_2` | GPL-3.0 | Chaumeil PA *et al.* 2022. https://doi.org/10.1093/bioinformatics/btac672 |
| geNomad | https://github.com/apcamargo/genomad | `genomad` 1.12.0, Bioconda build `pyhdfd78af_0` | BSD-4-Clause | Camargo AP *et al.* 2023. https://doi.org/10.1038/s41587-023-01953-y |
| CheckV | https://bitbucket.org/berkeleylab/checkv | `checkv` 1.1.1, Bioconda build `pyh106432d_0` | BSD-3-Clause-LBNL | Nayfach S *et al.* 2021. https://doi.org/10.1038/s41587-020-00774-7 |
| DIAMOND | https://github.com/bbuchfink/diamond | 2.1.11 in `checkm2`; 2.2.1 in `checkv_env` | GPL-3.0-or-later | Buchfink B *et al.* 2021. https://doi.org/10.1038/s41592-021-01101-x |
| Prodigal | https://github.com/hyattpd/Prodigal | 2.6.3; the exact build differs by locked environment | GPL-3.0-or-later | Hyatt D *et al.* 2010. https://doi.org/10.1186/1471-2105-11-119 |
| YaHS | https://github.com/c-zhou/yahs | `yahs` 1.2.2, Bioconda build `he4a0461_0` | MIT | Zhou C *et al.* 2023. https://doi.org/10.1093/bioinformatics/btac808 |
| bin3C Python 3 port | https://github.com/1001shiyuan/bin3C-python3 (fork of https://github.com/cerebis/bin3C) | Git commit `eb71dbf166e102e1c2fc2fa15f0112ad2fe60577` | AGPL-3.0 | DeMaere MZ, Darling AE. 2019. https://doi.org/10.1186/s13059-019-1643-1 |

## Notes on METAHICT-owned code and data

MetaCC and ImputeCC are maintained by the METAHICT authors and are not listed
as third-party dependencies here.  Their source remains attributed within the
repository.

Large reference databases are not distributed in the repository or container
image.  Users obtain them with the scripts in `installation/db/` and provide
their paths at run time.  The explicitly selected database releases are GTDB
release 220 and CheckV database v1.5.  The geNomad and CheckM2 download tools
record the database version in the downloaded database directory.

The bin3C fork retains upstream Git history and its AGPL-3.0 licence.  The
METAHICT environment installs it from the immutable commit above through
`installation/pip-requirements.txt`; no bin3C source tree is included in the
active METAHICT module directories.
