FROM community.wave.seqera.io/library/blast_diamond_hmmer_infernal_pruned:24ef3c0eea00bdb8

# Install build dependencies for C++ components (hyplas + unicycler)
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential python3-dev zlib1g-dev pigz && \
    rm -rf /var/lib/apt/lists/*

# Install read QC tools
RUN conda install -y -c bioconda -c conda-forge fastp chopper multiqc quast && conda clean -afy
RUN pip install --no-cache-dir NanoPlot

# Install hyplas dependencies
COPY external/unicycler-modified-for-hyplas /opt/unicycler-modified-for-hyplas
WORKDIR /opt/unicycler-modified-for-hyplas
RUN pip install --no-cache-dir .

# Build and install hyplas
COPY src /opt/hyplas/src
COPY Makefile configure /opt/hyplas/
WORKDIR /opt/hyplas
RUN ./configure --prefix=/usr/local && make -j$(nproc) && make install

RUN rm -rf /opt/unicycler-modified-for-hyplas /opt/hyplas
WORKDIR /
