# Snakemake workflow: Nanopore-amplicon

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.31.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/NP-metagenomics.svg?branch=master)](https://travis-ci.org/snakemake-workflows/NP-metagenomics)

This is a Snakemake workflow for the nanopore profiling of full-length 16S rRNA gene. Along with the workflow, one minimized version of compiled `guppy` and `Greengene 13.8` database are provided. The whole pipleline contains demultiplexing, taxonomy assignment and building the feature table.

## Usage

You can edit the path to `basecalled_dir` and `results_dir` in the `config.yaml`, and start snakemake workflow in the console.

    snakemake -p --use-conda -j 8 --resources mem_mb=50000

Or pass the arguments in the console.

    snakemake -p --use-conda -j 8 --resources mem_mb=50000 --config basecalled_dir=/path/to/basecalled_fqs results_dir=/path/to/results

we provide a wrapper function to start the live-mode of snakemake workflow. It also works in a finished nanopore run. In default setting, it will automatically exit if no new basecalled fastq files upload within 15 mins.

    ./live_snake -h

```
Note: The real-time mode of snakemake workflow in analyzing nanopore amplicon sequencing.

Usage: ./live_snake [-i -o -j -m -t -h]
  -i, --input     Required, the directory path for basecalled fastqs.
  -o, --output    Required, define the output directory path.
  -j, --jobs      Optional, define the thread number for analysis. Integer, default: 8.
  -m, --mem_mb    Optional, define the maximum memory usage(MB). Integer, default: 50000.
  -t, --wait      Optional, waiting seconds to exit. Default, 900, i.e. exit if no new fastq files upload within 15 mins.
  -h, --help      Optional, help message.

Example:
./live_snake -i /path/to/basecalled_fqs -o /path/to/results
./live_snake -i /path/to/basecalled_fqs -o /path/to/results -j 30 -m 200000 -t 1800

```

## Output

```
results
├── bioms_out
├── demultiplexed_fq
├── ONT-L7-GG.biom
└── ONT-L7-GG.txt
```

`bioms_out` directory holds the separate biom files for real-time batches of nanopore sequencers.  
`demultiplexed_fq` directory holds demultiplexed fastq files, which keep track of the barcode information. (It can be deleted when YOU MAKE SURE THERE NO NEW FASTQS.)  
`ONT-L7-GG.biom` and `ONT-L7-GG.txt` is the merged feature table in `biom` and `tsv` format, respectively.   

## Depolyment

### Step 1: Clone the workflow from the Github

    git clone https://github.com/yanhui09/NP-metagenomics.git

### Step 2: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake -np --use-conda -j 8 --resources mem_mb=50000 --config basecalled_dir=/path/to/basecalled_fqs results_dir=/path/to/results

Execute the workflow locally via

    snakemake -p --use-conda -j 8 --resources mem_mb=50000 --config basecalled_dir=/path/to/basecalled_fqs results_dir=/path/to/results


## Maintainence

### Step 1: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 2: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/NP-metagenomics.git` or `git remote add -f upstream https://github.com/snakemake-workflows/NP-metagenomics.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 3: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).