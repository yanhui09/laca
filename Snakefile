#--------------
configfile: "config.yaml"
#--------------
INPUT_DIR = config["basecalled_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
DATABASE_DIR = config["database_dir"].rstrip("/")

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all:
    input:
        OUTPUT_DIR + "/count_matrix.tsv",
        OUTPUT_DIR + "/sepp/ftable_filtered.biom",
        OUTPUT_DIR + "/sepp/sepp.tre",

include: "rules/demultiplex.smk"
#include: "rules/umi.smk"
include: "rules/umap_denoise.smk"
include: "rules/phylo.smk"