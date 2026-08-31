# METAHICT Conda locks

The Linux `.explicit.txt` files are the canonical release locks for the
environments used by METAHICT. Each line names an exact Conda package artifact
(channel, platform, version, and build), so `conda create --file` does not
re-solve dependency versions at installation time.

The human-readable YAML files in `installation/` remain the declared package
specifications. The explicit lock files are the reproducible installation
inputs for the tagged Linux release, including the separate MetaBAT2 runtime
that supplies `jgi_summarize_bam_contig_depths` without storing its binary in
the repository. `../pip-requirements.txt` separately pins Binning_refiner to a
checksum-verified upstream wheel and the METAHICT-compatible bin3C, ImputeCC,
and MetaCC packages to immutable Git commits.

`metahict_nextflow_env.explicit.txt` is the small Java environment that starts
Nextflow. It is intentionally separate from `metahict_env`, which contains the
scientific Python and bioinformatics runtime.

The installer hashes an existing Conda environment's explicit package list
before reusing it. If it does not match the release lock, installation stops
instead of silently running a different dependency set.

To refresh a lock for a deliberately tested environment, use the exact tested
environment and record the resulting SHA-256 values in `SHA256SUMS`:

```bash
conda list --explicit -p /path/to/environment > environment.explicit.txt
sha256sum environment.explicit.txt
```

Locks are platform-specific. A future macOS build would require separately
tested macOS locks; it could not reuse the Linux files.
