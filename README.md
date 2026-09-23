
# HyPlas
HyPlAs is a tool for assembling plasmids from hybrid short- and long-read sequencing data from bacterial isolates.
HyPlAs's main novelty is its use of a plasmid classification tool (such as platon) on short-read assembled contigs to guide plasmidic long-read selection and hybrid plasmid assembly.
HyPlAs is designed to work with single-genome sequencing data and has not been tested on metagenomic data.  

## Installation

HyPlAs is a C++ binary (`hyplas`) that orchestrates external bioinformatics
tools. Build it from source or use a container image.

### Build from source

Requires a C++20 compiler, CMake (>= 3.20), zlib, and network access at configure
time (CMake fetches [gtl](https://github.com/greg7mdp/gtl) and
[shrn](https://github.com/f0t1h/shrn). The runtime tools (Unicycler, Platon, minigraph,
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

HyPlAs is a pipeline that combines existing tools with the C++ `hyplas` binary.
It consists of the following steps (see the figure below): 
1. Read preprocessing;  
	1.a. Short reads are preprocessed with <a href="https://github.com/OpenGene/fastp">fastp</a>,  
    	1.b. Long reads are preprocessed with <a href="https://github.com/wdecoster/chopper">chopper</a>;  
2. Short reads are assembled using <a href="https://github.com/rrwick/Unicycler">Unicycler</a>;
3. Putative plasmidic long reads are detected in four stages:  
	3.a. The plasmid-contig classification tool <a href="https://github.com/oschwengers/platon">Platon</a> is used to detect plasmidic short-read contigs,  
	3.b. Long reads are mapped to the assembly graph using <a href="https://github.com/lh3/minigraph">minigraph</a>,  
   	3.c. Long-read mappings to short-read contigs and Platon results are used to select an initial set of putative plasmidic long reads,  
   	3.d. The set of putative plasmidic long reads is augmented by iteratively detecting overlapping long reads;  
4. The full short-read assembly graph generated in step 2 is refined with the plasmidic long reads selected during step 3, using Unicycler.

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
- `--platon-db`: Database used by Platon (<a href="https://zenodo.org/record/4066768/files/db.tar.gz">https://zenodo.org/record/4066768/files/db.tar.gz</a>)
```
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz
rm db.tar.gz
# Move the database to a suitable location
```
- `-s`: Space-separated short-read files
- `-l`: Long-read file (must be gzipped)
- `-o`: Output folder
- `-p`: Number of long-read recovery rounds to run (recommended: 2)
- `-t`: Number of threads (default: 16)
- `--sr-assembly`: Use a pre-computed short-read assembly graph instead of running Unicycler
- `--use-spades`: Assemble short reads with SPAdes instead of Unicycler
- `--per-component`: Bin contigs into graph components and run the hybrid assembly per component
- `--soft-fail`: On a failed hybrid assembly, fall back to the circular short-read contigs instead of exiting
- `--force`: Re-run every step (see below)
- `--keep-temp`: Keep each step's scratch directory under `<output>/tmp`

### Resuming a run

Re-running HyPlAs on an existing output folder redoes only work that is out of
date. Each step declares the files it reads and writes; a step is skipped when
all of its outputs exist and none is older than any of its inputs, just as
`make` decides. Changing an input (for example, re-running QC on the reads)
therefore re-runs the short-read assembly and everything downstream, while an
interrupted run picks up where it stopped. The covered steps are the short-read
assembly, Platon classification, minigraph alignment, initial read selection,
and missing-read extraction. Long-read recovery rounds and hybrid assemblies
are re-used whenever their output exists; pass `--force` to redo them, or
delete the corresponding `prop_lr/` or `unicycler_lr_*` entries.

`--force` re-runs every step regardless of timestamps. Scratch files for each
step live in a private directory under `<output>/tmp` and are removed when the
step finishes; `--keep-temp` keeps them, and the path is logged at the end of
each step that used one.

## Example run

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
HyPlAs creates the following files and directories in the output folder:  
- plasmids.final.it{iteration}.fasta:   
	- Assembled plasmids in FASTA format; iteration numbers 0 to {-p} are the results of each long-read recovery-round setting.
- unicycler_sr (directory):  
	- Short-read-only assembly by Unicycler;  
- classify (directory):  
	- classify/result.log: Platon log file,  
	- classify/result.json: Platon classification details for the short-read assembly contigs in JSON format,  
	- classify/result.tsv: Platon classification details for the short-read assembly contigs in TSV format,  
	- classify/result_p.tsv: List of contigs predicted by Platon as plasmidic or chromosomal;  
- lr2assembly.gaf: Graph alignment of long reads to the short-read-only assembly contigs;  
- plasmid_long_reads/plasmid.fastq.gz: Long reads classified as plasmidic by HyPlAs, in FASTQ format;  
- prop_lr/ (directory):
	- prop_lr/lr.round.[0-9]+.paf: Mappings of known plasmid long reads to unknown long reads; the integer suffix indicates the iteration round of plasmidic long-read augmentation (step 3.d),  
	- prop_lr/lr.round.[0-9]+.fastq.gz: Plasmidic long-read sequences recovered in augmentation round X (X in [0-9]);  
- unicycler_lr_{iteration} (directories):  
	- Outputs of the final Unicycler hybrid assembly using the short-read-only assembly and predicted plasmidic long reads for each iteration.  
