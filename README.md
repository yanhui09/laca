# Kamp

[![Snakemake](https://img.shields.io/badge/snakemake-=7.8.5-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/Kamp.svg?branch=master)](https://travis-ci.org/snakemake-workflows/Kamp)

`Kamp` is a kmer-based denoise workflow to analyze long-read noisy amplicons by Nanopore sequencing, e.g., 16S rRNA gene.
Using `Snakemake` as a job controller behind, `Kamp` is wrapped into a python package for development and mainteniance.
`Kamp` provides an end-to-end solution from bascecalled reads to the final count matrix.

**Important: `Kamp` is under development, and here released as a preview for early access. 
As `Kamp` integrates multiple bioinformatic tools, it is only tested in Linux systems, i.e., Ubuntu.
The detailed instruction and manuscript is in preparation.**

# Installation
[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) is the only required dependency prior to installation.
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) is enough for the whole pipeline. 

1. Clone the Github repository and create an isolated `conda` environment
```
git clone https://github.com/yanhui09/Kamp.git
cd Kamp
conda env create -n kamp -f kampenv.yml 
```
You can speed up the whole process if [`mamba`](https://github.com/mamba-org/mamba) is installed.
```
mamba env create -n kamp -f kampenv.yml 
```
2. Install `Kamp` with `pip`
We suggest installing `Kamp` in the above `conda` environment
```
conda activate kamp
pip install --editable .
```

At this moment, `Kamp` uses a compiled but tailored `guppy` for demultiplexing the barcodes sets in our lab.
**Remember to revise the config file in `guppy` if new barcodes are introduced.**

# Example
```
conda activate kamp                                  # activate required environment 
kamp init -f /path/to/fastqs -d /path/to/database    # init config file
kamp run all                                         # start analysis
```

# Usage

```
Usage: kamp [OPTIONS] COMMAND [ARGS]...

  Kamp: a k-mer based denoise pipeline to process long read amplicon
  sequencing by Nanopore. To follow updates and report issues, see:
  https://github.com/yanhui09/Kamp.

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  init  Prepare the config file.
  run   Run Kamp workflow.
```

`Kamp` is easy to use. You can start a new analysis in two steps using `kamp init` and `kamp run` . 

Remember to first activate the conda environment.
```
conda activate kamp
```

1. Intialize a config file with `kamp init`

`kamp init` will generate a config file in the working directory, which contains the necessary parameters to run `Kamp`.

```
Usage: kamp init [OPTIONS]

  Prepare the config file for Kamp.

Options:
  -f, --fqdir PATH                Path to the basecalled fastq files.
                                  [required]
  -d, --dbdir PATH                Path to the taxonomy databases.  [required]
  -w, --workdir PATH              Output directory for Kamp.  [default: .]
  --fqs-min INTEGER               Minimum number of reads for the
                                  demultiplexed fastqs.  [default: 1000]
  --no-pool                       Do not pool the reads for denoising.
  --subsample                     Subsample the reads.
  --no-trim                       Do not trim the primers.
  --kmerbin                       Conduct kmer binning.
  --cluster [kmerCon|clustCon|isONclustCon|isONcorCon|umiCon]
                                  Methods to generate consensus.  [default:
                                  isONclustCon]
  --chimerf                       Filter possible chimeras by vsearch.
  --jobs-min INTEGER              Number of jobs for common tasks.  [default:
                                  2]
  --jobs-max INTEGER              Number of jobs for threads-dependent tasks.
                                  [default: 6]
  -h, --help                      Show this message and exit.
```

2. Start analysis with `kamp run`

`kamp run` will trigger the full workflow or a specfic module under defined resource accordingly.
Get a dry-run overview with `-n`. `Snakemake` arguments can be appened to `kamp run` as well.

```
Usage: kamp run [OPTIONS] {demultiplex|qc|kmerBin|kmerCon|clustCon|isONclustCo
                n|isONcorCon|umiCon|quant|taxa|tree|requant|all|initDB|nanosim
                } [SNAKE_ARGS]...

  Run Kamp workflow.

Options:
  -w, --workdir PATH  Output directory for Kamp.  [default: .]
  -j, --jobs INTEGER  Maximum jobs to run in parallel.  [default: 6]
  -m, --maxmem FLOAT  Maximum memory to use in GB.  [default: 50]
  -n, --dryrun        Dry run.
  -h, --help          Show this message and exit.
```
