# LACA: Long Amplicon Consensus Analysis

[![snakemake](https://img.shields.io/badge/snakemake-7.22.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![linux/amd64](https://img.shields.io/badge/linux-amd64-blue.svg)](https://en.wikipedia.org/wiki/X86-64)
[![CI](https://github.com/yanhui09/laca/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/yanhui09/laca/actions?query=branch%3Amaster+workflow%3ACI)
[![Docker](https://github.com/yanhui09/laca/actions/workflows/docker.yml/badge.svg?branch=master)](https://github.com/yanhui09/laca/actions?query=branch%3Amaster+workflow%3ADocker)

`LACA` is a reproducible and scalable workflow for **Long Amplicon Consensus Analysis**, e.g., 16S rRNA gene.
Using `snakemake` as the job controller, `LACA` is wrapped into a python package for development and mainteniance.
`LACA` provides an end-to-end solution from bascecalled reads to the final count matrix.

**Important: 
As `LACA` integrates multiple bioinformatic tools, it is only tested in Linux systems, i.e., Ubuntu.
**

[preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2023.11.26.568539v1)

# Docker image
The easiest way to use `LACA` is to pull the `docker` image from [Docker Hub](https://hub.docker.com/r/yanhui09/laca) for cross-platform support.
```
docker pull yanhui09/laca
```

**To use the docker image**, you need to mount your data directory, e.g., `pwd`, to the  `/home` in the container.
```
docker run -it -v `pwd`:/home --privileged yanhui09/laca
```

# Installation from GitHub repository
[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) is the only required dependency prior to installation.
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) is enough for the whole pipeline. 

1. Clone the Github repository and create an isolated `conda` environment
```
git clone https://github.com/yanhui09/laca.git
cd laca
conda env create -n laca -f env.yaml 
```
You can speed up the whole process if [`mamba`](https://github.com/mamba-org/mamba) is installed.
```
mamba env create -n laca -f env.yaml 
```
2. Install `LACA` with `pip`
      
To avoid inconsistency, we suggest installing `LACA` in the above `conda` environment
```
conda activate laca
pip install --editable .
```

At this moment, `LACA` uses a compiled but tailored `guppy` for barcode demultiplexing (in our lab).<br>
**Remember to prepare the barcoding files in `guppy` if new barcodes are introduced.** [Click me](laca/workflow/resources/README.md)

# Example
```
laca init -b /path/to/basecalled_fastqs -d /path/to/database    # init config file and check
laca run all                                         # start analysis
```

# Usage

```
Usage: laca [OPTIONS] COMMAND [ARGS]...

  LACA: a reproducible and scaleable workflow for Long Amplicon Consensus
  Analysis. To follow updates and report issues, see:
  https://github.com/yanhui09/laca.

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  init  Prepare the config file.
  run   Run LACA workflow.
```

`LACA` is easy to use. You can start a new analysis in two steps using `laca init` and `laca run` . 

Remember to activate the `conda` environment if `LACA` is installed in a `conda` environment.
```
conda activate laca
```

1. Intialize a config file with `laca init`

`laca init` will generate a config file in the working directory, which contains the necessary parameters to run `LACA`.

```
Usage: laca init [OPTIONS]

  Prepare config file for LACA.

Options:
  -b, --bascdir PATH              Path to a directory of the basecalled fastq
                                  files. Option is mutually exclusive with
                                  'merge', 'merge_parent', 'demuxdir'.
  -x, --demuxdir PATH             Path to a directory of demultiplexed fastq
                                  files. Option is mutually exclusive with
                                  'merge', 'merge_parent', 'bascdir'.
  --merge PATH                    Path to the working directory of a completed
                                  LACA run  [Mutiple]. Runs will be combined
                                  if --merge_parent applied. Option is
                                  mutually exclusive with 'bascdir',
                                  'demuxdir'.
  --merge-parent PATH             Path to the parent of the working
                                  directories of completed LACA runs. Runs
                                  will be combined if --merge applied. Option
                                  is mutually exclusive with 'bascdir',
                                  'demuxdir'.
  -d, --dbdir PATH                Path to the taxonomy databases.  [required]
  -w, --workdir PATH              Output directory for LACA.  [default: .]
  --demuxer [guppy|minibar]       Demultiplexer.  [default: guppy]
  --fqs-min INTEGER               Minimum number of reads for the
                                  demultiplexed fastqs.  [default: 1000]
  --no-pool                       Do not pool the reads for denoising.
  --subsample                     Subsample the reads.
  --no-chimera-filt               Do not filter chimeric reads.
  --no-primer-check               Do not check primer pattern.
  --cluster [isONclust|umapclust|meshclust]
                                  Cluster approaches.  [Mutiple]  [default:
                                  isONclust, meshclust]
  --consensus [kmerCon|miniCon|isoCon|umiCon]
                                  Consensus methods.  [Mutiple]  [default:
                                  kmerCon]
  --quant [seqid|minimap2]        Create abundance matrix by sequence id or
                                  minimap2.  [Mutiple]  [default: seqid]
  --uchime                        Filter chimeras by uchime-denovo in vsearch.
  --jobs-min INTEGER              Number of jobs for common tasks.  [default:
                                  2]
  --jobs-max INTEGER              Number of jobs for threads-dependent tasks.
                                  [default: 6]
  --ont                           Use config template for ONT reads. Option is
                                  mutually exclusive with 'isoseq'.
  --isoseq                        Use config template for PacBio CCS reads.
                                  Option is mutually exclusive with 'ont'.
  --longumi                       Use primer design from longumi paper (https:
                                  //doi.org/10.1038/s41592-020-01041-y).
                                  Option is mutually exclusive with
                                  'simulate'.
  --simulate                      Use config template for in silicon test.
                                  Option is mutually exclusive with 'longumi'.
  --clean-flags                   Clean flag files.
  -h, --help                      Show this message and exit.
```

2. Start analysis with `laca run`

`laca run` will trigger the full workflow or a specfic module under defined resource accordingly.
Get a dry-run overview with `-n`. `Snakemake` arguments can be appened to `laca run` as well.

```
Usage: laca run [OPTIONS] {demux|qc|clust|kmerCon|miniCon|isoCon|umiCon|quant|
                taxa|tree|all|merge|initDB|simulate} [SNAKE_ARGS]...

  Run LACA workflow.

Options:
  -w, --workdir PATH     Working directory for LACA.  [default: .]
  -c, --configfile FILE  Config file for LACA. Use config.yaml in working
                         directory if not specified.
  -j, --jobs INTEGER     Maximum jobs to run in parallel.  [default: 6]
  -m, --maxmem FLOAT     Specify maximum memory (GB) to use. Memory is
                         controlled by profile in cluster execution.
  --profile TEXT         Snakemake profile for cluster execution.
  -n, --dryrun           Dry run.
  -h, --help             Show this message and exit.
```