# Snakemake workflow: Kamp

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.13.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/Kamp.svg?branch=master)](https://travis-ci.org/snakemake-workflows/Kamp)

This is a kmer-based denoise pipeline for the nanopore amplicon sequencing, e.g., 16S rRNA gene. Along with the workflow, one minimized version of compiled `guppy` was provided for demultiplexing.

# Installation
Make sure conda is installed. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is enough for the whole pipeline.

Install latest version of mamba and snakemake, and activate the conda environment to run snakemake.
```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

# Usage

You can edit the path to `basecalled_dir` and `results_dir` in the `config.yaml`, and start snakemake workflow in the console.
```
snakemake -p --use-conda -j 8 --resources mem_mb=50000
```

Or pass the arguments in the console.
```
snakemake -p --use-conda -j 8 --resources mem_mb=50000 --config basecalled_dir=/path/to/basecalled_fqs results_dir=/path/to/results threads_taxa=8
```

Initialize the database setup
```
snakemake --core 6 --use-conda init_database
```

