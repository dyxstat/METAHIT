# Publishing the METAHICT release container

The released workflow will use one all-tools image.  It contains the exact
locked Conda environments `metahict_env`, `checkm2`, `gtdbtk-2.4.0`,
`genomad`, and `ccfind_env`.  Large reference databases remain external user-configured
paths and are not included in the image.

The planned registry location is:

```text
ghcr.io/1001shiyuan/metahict:<tag-you-choose-at-publication>
```

This is only a planned location until the image is pushed and made public.

## Release sequence

From the METAHICT repository root:

```bash
# Build locally; this does not upload anything.
nextflow/containers/release_image.sh build <test-tag>
nextflow/containers/release_image.sh inspect <test-tag>
```

Run the complete `example_data` workflow using the local image and the Docker
profile.  The single-image override ensures that every process uses exactly
the candidate release image:

```bash
nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile docker \
  --container_image_override ghcr.io/1001shiyuan/metahict:<test-tag> \
  --samplesheet nextflow/assets/example_data_samplesheet.csv
```

After all Docker and Apptainer tests are successful, authenticate and publish:

```bash
docker login ghcr.io
nextflow/containers/release_image.sh push <tag-you-choose-at-publication>
```

The push command prints an immutable image digest.  Record that digest in the
tagged release configuration, for example:

```text
--container_image_override \
  ghcr.io/1001shiyuan/metahict@sha256:<published-digest>
```

Finally, set the GitHub Container Registry package visibility to **Public**.
Verify that a clean machine can pull the digest without logging in, then run
the complete minimal example under both Docker and Apptainer.  Record the
tested tag and digest in the release notes and permanent archive.
