# METAHI-T Nextflow Release Checklist

This checklist is for releasing the native DSL2 workflow layer. It does not
change METAHI-T module behavior.

## Release identity

1. Add the completed changes to `nextflow/CHANGELOG.md`.
2. Record the immutable container image digest or SIF filename used for the release.
3. Record the METAHI-T commit or archive that the workflow was tested against.
4. Create an annotated Git tag on that exact commit.  The release identifier is
   chosen by the maintainers at release time; do not create or change one while
   work is still in progress.

## Tests

1. Run the local orchestration smoke test:

   ```bash
   bash nextflow/ci/run_smoke_ci.sh
   ```

2. Run the containerized example-data workflow on a server with databases:

   ```bash
   nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
     -profile apptainer \
     --samplesheet nextflow/assets/example_data_samplesheet.csv \
     --out_root /raid/projects/Shiyuan/metahict/nextflow_results \
     --report_dir /raid/projects/Shiyuan/metahict/nextflow/reports/example_data_apptainer \
     -work-dir /raid/projects/Shiyuan/metahict/nextflow/work/example_data_apptainer \
     --threads 80 \
     --clean true \
     -ansi-log false
   ```

3. Validate expected outputs:

   ```bash
   python nextflow/bin/check_expected_outputs.py \
     --root /raid/projects/Shiyuan/metahict/nextflow_results \
     --manifest nextflow/tests/expected/example_data_outputs.tsv
   ```

## Release Notes

The release notes should state:

- Nextflow version used for testing.
- Container engine used for testing.
- Immutable container image digest or SIF filename.
- Example dataset size.
- Which modules were tested end to end.

## Permanent archive

After all tests pass and the annotated Git tag has been created:

1. Create a GitHub Release from that exact tag.
2. Archive that GitHub Release with an archival service such as Zenodo.
3. Record the resulting permanent DOI, Git tag, commit SHA, and container
   image digest in the GitHub Release notes, manuscript Code Availability
   statement, and Supplementary Information.
4. Confirm that the archived record links to the public source repository and
   that the container image can be pulled without authentication.

No DOI, Git tag, or release identifier is recorded in this working tree before
the final tested release exists.
