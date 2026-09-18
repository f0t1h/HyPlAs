FROM community.wave.seqera.io/library/blast_diamond_hmmer_infernal_pruned:24ef3c0eea00bdb8

# Install build dependencies for C++ components (hyplas + unicycler)
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential cmake git python3-dev zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Install tools used by the core Nextflow workflow. QUAST remains in the
# separate post-run analysis environment.
RUN micromamba install -y -n base --root-prefix /opt/conda --strict-channel-priority \
        -c conda-forge -c bioconda \
        "python>=3.10,<3.13" pip biopython spades minigraph minimap2 racon \
        platon prodigal mummer4 fastp chopper pigz multiqc && \
    micromamba clean -afy

# Install hyplas dependencies
COPY external/unicycler-modified-for-hyplas /opt/unicycler-modified-for-hyplas
WORKDIR /opt/unicycler-modified-for-hyplas
RUN pip install --no-cache-dir .

# Build and install hyplas
COPY src /opt/hyplas/src
COPY CMakeLists.txt /opt/hyplas/
WORKDIR /opt/hyplas
RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
        -DBUILD_TESTING=OFF && \
    cmake --build build --parallel 2 && \
    cmake --install build && \
    test -x /usr/local/bin/hyplas

RUN hyplas --check-deps
RUN command -v fastp && command -v chopper && command -v multiqc

RUN rm -rf /opt/unicycler-modified-for-hyplas /opt/hyplas
WORKDIR /
