# METAHICT Conda locks

The `linux-64/*.explicit.txt` files are the canonical release locks for the
environments used by METAHICT. Each line names an exact Conda package artifact
(channel, platform, version, and build), so `conda create --file` does not
re-solve dependency versions at installation time.

The human-readable YAML files in `installation/` remain the declared package
specifications. The explicit lock files are the reproducible installation
inputs for the tagged Linux release. `../pip-requirements.txt` separately pins
the non-Conda METAHICT-compatible bin3C, ImputeCC, and MetaCC packages at
immutable Git commits.

The installer hashes an existing Conda environment's explicit package list
before reusing it. If it does not match the release lock, installation stops
instead of silently running a different dependency set.

To refresh a lock for a deliberately tested environment, use the exact tested
environment and record the resulting SHA-256 values in `SHA256SUMS`:

```bash
conda list --explicit -p /path/to/environment > environment.explicit.txt
sha256sum environment.explicit.txt
```

Locks are platform-specific. A macOS build must be tested and released with
separate `osx-*` locks; it must not reuse the Linux files.
