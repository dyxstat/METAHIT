FROM condaforge/miniforge3:24.11.3-0

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG METAHICT_ENVS="metahict_env"
ARG METAHICT_VERSION="unreleased"
ARG VCS_REF="unknown"
ARG SOURCE_REPOSITORY="https://github.com/dyxstat/METAHICT"
LABEL org.opencontainers.image.title="METAHICT" \
      org.opencontainers.image.description="METAHICT native DSL2 workflow runtime" \
      org.opencontainers.image.source="${SOURCE_REPOSITORY}" \
      org.opencontainers.image.version="${METAHICT_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.licenses="GPL-3.0-only"

ENV CONDA_BASE=/opt/conda
ENV PATH=/opt/conda/envs/metahict_env/bin:/opt/conda/bin:/usr/local/bin:${PATH}

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        bzip2 \
        coreutils \
        curl \
        findutils \
        gawk \
        git \
        grep \
        gzip \
        procps \
        sed \
        tar \
        unzip && \
    rm -rf /var/lib/apt/lists/*

COPY installation/locks/linux-64 /tmp/metahict-locks
COPY installation/pip-requirements.txt /tmp/pip-requirements.txt

RUN conda config --set channel_priority flexible && \
    for env_name in ${METAHICT_ENVS}; do \
        case "${env_name}" in \
            metahict_env) \
                conda create -y -n metahict_env --file /tmp/metahict-locks/metahict_env.explicit.txt ;; \
            checkm2) \
                conda create -y -n checkm2 --file /tmp/metahict-locks/checkm2.explicit.txt ;; \
            gtdbtk-2.4.0) \
                conda create -y -n gtdbtk-2.4.0 --file /tmp/metahict-locks/gtdbtk-2.4.0.explicit.txt ;; \
            genomad) \
                conda create -y -n genomad --file /tmp/metahict-locks/genomad.explicit.txt ;; \
            ccfind_env) \
                conda create -y -n ccfind_env --file /tmp/metahict-locks/ccfind_env.explicit.txt ;; \
            *) \
                echo "Unknown METAHICT env requested: ${env_name}" >&2; exit 1 ;; \
        esac; \
    done && \
    conda clean -afy

RUN if [ -x /opt/conda/envs/metahict_env/bin/python ]; then \
        /opt/conda/envs/metahict_env/bin/python -m pip install --no-deps --no-cache-dir \
            -r /tmp/pip-requirements.txt && \
        /opt/conda/envs/metahict_env/bin/python -c "import bin3c_python3"; \
    fi

RUN printf '#!/usr/bin/env bash\nexit 0\n' > /usr/local/bin/module && \
    chmod +x /usr/local/bin/module && \
    for spec in \
        checkm2:checkm2 \
        gtdbtk:gtdbtk-2.4.0 \
        genomad:genomad; do \
        tool="${spec%%:*}"; env_name="${spec##*:}"; \
        printf '#!/usr/bin/env bash\nexec conda run -n %s %s "$@"\n' "${env_name}" "${tool}" > "/usr/local/bin/${tool}"; \
        chmod +x "/usr/local/bin/${tool}"; \
    done

WORKDIR /workspace
