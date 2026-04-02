FROM community.wave.seqera.io/library/blast_diamond_hmmer_infernal_pruned:24ef3c0eea00bdb8

# Install build dependencies for unicycler C++ components
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential python3-dev zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy and install the modified unicycler
COPY external/unicycler-modified-for-hyplas /opt/unicycler-modified-for-hyplas
WORKDIR /opt/unicycler-modified-for-hyplas
RUN pip install --no-cache-dir .

WORKDIR /
