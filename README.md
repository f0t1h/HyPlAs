
# HyPlas
HyplAs is a tool aimed at assembling plasmids from hybid short-read and long-read sequencing data for bacerial isolates.
HyPlAs main novlty is to incorporate a plasmid classification tools (Such as platon) on short-read assembled contigs to aid plasmidic long-read selection and performs hybrid plasmids assembly.
HyPlAs has been desiged to work with single genome sequencing data, and has not been tested on metagenomics data.  

## Installation

HyPlAs is a C++ binary (`hyplas`) that orchestrates external bioinformatics
tools. Build it from source or use a container image.

### Build from source

Requires a C++20 compiler, CMake (>= 3.20), zlib, and network access at configure
time (CMake fetches [gtl](https://github.com/greg7mdp/gtl) and
[shrn](https://github.com/f0t1h/shrn); a sibling `../shrn` checkout is used
automatically if present). The runtime tools (Unicycler, Platon, minigraph,
minimap2, SPAdes, ...) are listed in `environment.yml` and can be installed
with conda/mamba:

```bash
mamba env create -f environment.yml
conda activate hyplas-env

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build            # installs into the active env prefix by default
ctest --test-dir build           # optional: run the unit tests
hyplas --help
hyplas --check-deps              # verify all external tools are on PATH
```

The container and `environment.yml` include everything used by `main.nf`:
Fastp, Chopper, HyPlAs, and MultiQC. QUAST is optional post-run analysis and
is kept in a separate environment:

```bash
mamba env create -f environment.analysis.yml
conda activate hyplas-analysis
```

### Container images

HyPlAs containers are built on top of pinned Seqera Wave base images that bundle the bioinformatics dependencies (BLAST, Diamond, HMMER, Infernal, etc.).

Build locally:
```bash
docker build -t hyplas:latest .
apptainer build hyplas.sif Apptainer.def
```

## Overview

HyPlAs is a pipeline combining existing tools and the C++ `hyplas` binary.
HyPlAs is composed of
the following steps (see figure below): 
1. Reads preprocessing;  
	1.a. Short reads are preprocessed with <a href="https://github.com/OpenGene/fastp">fastp</a>,  
    	1.b. Long reads are preprocessed with <a href="https://github.com/wdecoster/chopper">chopper</a>;  
3. Short reads are assembled using <a href="https://github.com/rrwick/Unicycler">Unicycler</a>;
4. The detection of putative plasmidic long reads is done in four stages:  
	3.a. the plasmid contigs classification tool <a href="https://github.com/oschwengers/platon">Platon</a> is used to detect plasmidic short-read contigs,  
	3.b. long reads are mapped to the assembly graph using <a href="https://github.com/lh3/minigraph">minigraph</a>,  
   	3.c. long-read mapping to short-read contigs and platon results are used to select an initial set of putative plasmidic long reads,  
   	3.d. the set of putative plasmidic long reads is augmented by iteratively detecting overlapping long reads;  
5. The full short-read assembly graph generated in step 2 is refined with the plasmidic long reads selected during step 3, using Unicycler.

![HyPlAs](resources/HyPlAs_pipeline.png?raw=true)

## Usage

Run the pipeline directly:
```bash
hyplas --platon-db db -s sr_1.fastq sr_2.fastq -l lr.fastq.gz -o hyplas-out -t 16 -p 2
```

Or run the full QC + assembly + reporting workflow with Nextflow (see `main.nf`):
```bash
nextflow run main.nf --samples samplesheet.csv --platonDb db --outdir results
```
### Input
 - --platon-db: Database used by Platon (<a href="https://zenodo.org/record/4066768/files/db.tar.gz">https://zenodo.org/record/4066768/files/db.tar.gz</a>)
 ```
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz
rm db.tar.gz
# Move the database to a suitable location
```
- -s space separated short read files
- -l long reads file (required to be gzipped)
- -o output folder
- -p number of long-read recovery rounds to be executed (Recommend 2 rounds)

## Example run.

### Download long and short reads from SRA (SAMN05238672)
```
fasterq-dump SRR3666207 SRR10173103
```

### Run QC tools
```
fastp --in1 SRR3666207_1.fastq --in2 SRR3666207_2.fastq --out1 SRR3666207_1.qc.fastq --out2 SRR3666207_2.qc.fastq --unpaired1 SRR3666207_unpaired.qc.fastq
chopper  -q 9 -l 500 --headcrop 75 --tailcrop 75 --input SRR10173103.fastq --threads 16 | gzip > SRR10173103.qc.fastq.gz 
```

### Download platon database
```
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz
rm db.tar.gz
```

### Run HyPlAs
```
hyplas -l SRR10173103.qc.fastq.gz -s SRR3666207_1.qc.fastq SRR3666207_2.qc.fastq -p 2 -o hyplas_outdir --platon-db db -t 64
```

### Output
HyPlAs creates in the output folder the following files and directories:  
- plasmids.final.it{iteration}.fasta:   
	- assembled plasmids, in FASTA format; iteration numbers 0 to {-p} results of each long-read recovery rounds settings.
- unicycler_sr (directory):  
	- short-read-only assembly by Unicycler;  
- classify (directory):  
	- classify/result.log: Platon log file,  
	- classify/result.json: Platon classification details of the short-read assembly contigs in json format,  
	- classify/result.tsv: Platon classification details of the short-read assembly contigs in tsv format,  
	- classify/result_p.tsv: List of contigs predicted by Platon as plasmidic or chromosomal;  
- lr2assembly.gaf: Graph alignment of long reads to the short-read-only assembly cotigs;  
- plasmid_long_reads/plasmid.fastq.gz: Long-reads classified as plasmidic by HyPlAs, in FASTQ format;  
- prop_lr/ (directory):
	- prop_lr/lr.round.[0-9]+.paf: Mappings of the known plasmid long-reads to unknown long-reads, the integer suffix indicates the iteration round of plasmidic long-reads augmentation (step 3.d),  
	- prop_lr/lr.round.[0-9]+.fastq.gz: Plasmidic long-read sequences recovered in augmentation round X (X in [0-9]);  
- unicycler_lr_{iteration} (directories):  
	- Outputs of the final Unicycler hybrid assembly using the short-read--only assembly and the predicted plasmdic long reads for each iteration.  
