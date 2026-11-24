FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba create -y -p /opt/conda/envs/bio \
      -c conda-forge -c bioconda \
      python=3.10 pip tini openjdk=17 bash gzip \
      bwa=0.7.19 samtools=1.25 minimap2=2.31 \
      trimmomatic=0.39 picard=3.5.0 msamtools=1.1.3 \
      prinseq-plus-plus=1.2.4 nanofilt=2.8.0 \
    && micromamba create -y -p /opt/conda/envs/app \
      -c conda-forge python=3.11 pip \
    && micromamba clean -afy

COPY mapping_tracker/ /app/mapping_tracker/
COPY mapped_reference_extract/extract_remap.py /app/mapped_reference_extract/extract_remap.py

WORKDIR /app/mapping_tracker

RUN /opt/conda/envs/app/bin/python -m pip install --no-cache-dir -r requirements.txt \
    && mkdir -p /home/mapping \
      /app/mapping_tracker/app/cache/browser \
      /app/mapping_tracker/mapping_results \
      /app/mapping_tracker/references \
      /app/mapping_tracker/samples \
      /data \
    && chown -R 1000:1000 /home/mapping \
      /app/mapping_tracker/app/cache \
      /app/mapping_tracker/mapping_results \
      /app/mapping_tracker/references \
      /app/mapping_tracker/samples \
      /data

ENV PATH=/opt/conda/envs/app/bin:/opt/conda/envs/bio/bin:$PATH
ENV JAVA_HOME=/opt/conda/envs/bio
ENV HOME=/home/mapping
ENV MAPPING_DEDUP_HEAP=8g

USER 1000:1000

ENTRYPOINT ["/opt/conda/envs/bio/bin/tini", "--"]
CMD ["streamlit", "run", "app/app.py", "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
