# METAHICT Nextflow Containers

This folder adds optional container support for the native DSL2 workflow. It
does not replace or modify the METAHICT module scripts.

The design is module-level containerization. Each image contains the named
conda environments needed by that module. For example, the Module 7 image
contains `metahict_env` and `checkm2`, so the existing reassembly script can
still switch between those environments inside the container.

Databases are not baked into the images. Large databases such as GTDB-Tk,
geNomad, CheckV, and CheckM2 should stay on shared storage and be supplied with
the usual METAHICT database path parameters.

## Build Docker Images

Run from the METAHICT repository root:

```bash
nextflow/containers/build_containers.sh docker
```

Build only one module image:

```bash
nextflow/containers/build_containers.sh docker binning
```

## Build Apptainer/Singularity Images

Run from the METAHICT repository root:

```bash
nextflow/containers/build_containers.sh apptainer
```

Images are written to:

```text
nextflow/containers/images/
```

## Run With Containers

Docker:

```bash
nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf -profile docker ...
```

Apptainer:

```bash
nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf -profile apptainer ...
```

The default `local` profile does not use containers and preserves the existing
repo-local conda environment behavior.

## Release Image

The release uses one image containing all five locked environments, rather
than publishing a separate image for each pipeline stage.  This is simpler for
users and still lets each DSL2 process select its own declared resource
requirements.  The eventual public location is `ghcr.io/1001shiyuan/metahict`.

Build and inspect the candidate image without uploading it:

```bash
nextflow/containers/release_image.sh build <test-tag>
nextflow/containers/release_image.sh inspect <test-tag>
```

See [PUBLISHING.md](PUBLISHING.md) for the required test, publication, digest,
and public-access checks.  Do not substitute a mutable tag for the final
digest in a tagged METAHICT release.
