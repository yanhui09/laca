# Snakemake workflow: Nanopore-amplicon

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.31.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/NP-metagenomics.svg?branch=master)](https://travis-ci.org/snakemake-workflows/NP-metagenomics)

This is a Snakemake workflow for the nanopore profiling of full-length 16S rRNA gene. Along with the workflow, one minimized version of compiled `guppy` and `Greengene 13.8` database are provided. The whole pipleline contains demultiplexing, taxonomy assignment and building the feature table.

## Usage

You can edit the path to `basecalled_dir` and `results_dir` in the `config.yaml`, and start snakemake workflow in the console.

    snakemake -p --use-conda -j 8 --resources mem_mb=50000

Or pass the arguments in the console.

    snakemake -p --use-conda -j 8 --resources mem_mb=50000 --config basecalled_dir=/path/to/basecalled_fqs results_dir=/path/to/results threads_taxa=8

