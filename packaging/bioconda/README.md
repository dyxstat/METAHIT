# Bioconda release recipe

This directory contains a locally buildable recipe for the METAHICT workflow
layer. It installs the `metahict` and `metahict-nextflow` commands and declares
Nextflow as a normal package dependency. The recipe test runs the installed
13-process DSL2 graph in stub mode and checks its expected output contract;
the real-data release test remains a separate prerequisite.

Before opening a pull request to `bioconda-recipes`, release maintainers must:

1. create and test an immutable METAHICT tag;
2. replace the recipe's local `source.path` with that tag's archive URL and
   SHA-256;
3. document that this is a thin workflow package and that the supported full
   runtime is installed from the repository's exact Conda locks;
4. run `conda-build`/Bioconda lint and the bundled example-data workflow; and
5. submit the recipe to Bioconda and use the package name only after it is
   accepted.

The repository-local exact Conda locks are the supported complete runtime
installation. This staging recipe must not be described as an already
published Bioconda package.
