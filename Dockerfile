FROM mambaorg/micromamba:1.5.8-jammy

LABEL org.opencontainers.image.title="RadarMeta"
LABEL org.opencontainers.image.description="Reconciled, ML-augmented metagenomics classification"
LABEL org.opencontainers.image.source="https://github.com/Ifeanyi-omah/radarmeta"
LABEL org.opencontainers.image.licenses="MIT"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps wget ca-certificates curl gawk bash gzip tar \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"

WORKDIR /pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER . /pipeline/
RUN chmod +x /pipeline/radarmeta /pipeline/bin/*.py

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
