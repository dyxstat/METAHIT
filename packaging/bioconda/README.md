# Bioconda release recipe

METAHICT is distributed through Bioconda. This directory contains the
repository-local staging copy of its recipe; the authoritative published
recipe is maintained in the `bioconda/bioconda-recipes` repository. It installs
the `metahict` and `metahict-nextflow` commands and declares Nextflow as a
normal package dependency. The recipe test runs the installed DSL2 graph in
stub mode and checks its expected output contract; the real-data release test
remains a separate prerequisite.

Before opening a pull request to `bioconda-recipes`, release maintainers must:

1. create and test an immutable METAHICT tag;
2. update the Bioconda recipe to the tagged archive and its SHA-256;
3. document that this is a thin workflow package and that the supported full
   runtime is installed from the repository's exact Conda locks;
4. run `conda-build`/Bioconda lint and the bundled example-data workflow; and
5. submit the recipe update to `bioconda/bioconda-recipes`.

The repository-local exact Conda locks are the supported complete runtime
installation.
