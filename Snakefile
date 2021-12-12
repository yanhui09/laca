#--------------
configfile: "config.yaml"
#--------------
INPUT_DIR = config["basecalled_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
DATABASE_DIR = config["database_dir"].rstrip("/")

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

def check_val_pool(val):    
    if not isinstance(val, bool):
        raise ValueError('\n\tPooling only allows bool type [True/False].\n\t{} is used in the config file.'.format(val))

rule all:
    input:
        OUTPUT_DIR + "/demultiplex_DONE",
        OUTPUT_DIR + "/qc_DONE",
        OUTPUT_DIR + "/kmerClust_DONE",
        OUTPUT_DIR + "/isONclustConsensus_DONE",
        OUTPUT_DIR + "/quant_DONE",
        OUTPUT_DIR + "/phylo_DONE",

include: "rules/demultiplex.smk"
#include: "rules/umi.smk"
include: "rules/qc.smk"
include: "rules/umap_hdbscan.smk"
include: "rules/isONclustCon.smk"
#include: "rules/isONcorrect_IsoCon.smk"
include: "rules/quant.smk"
include: "rules/phylo.smk"

rule demultiplex:
    input: lambda wc: expand(OUTPUT_DIR + "/raw/{barcode}.fastq", barcode=get_demultiplexed(wc))
    output: temp(touch(OUTPUT_DIR + "/demultiplex_DONE"))

rule qc:
    input: lambda wc: get_filt(wc, pooling = config["pooling"]),
    output: temp(touch(OUTPUT_DIR + "/qc_DONE"))

rule kmerClust:
    input: lambda wc: get_kmerClust(wc, pooling = config["pooling"]),
    output: temp(touch(OUTPUT_DIR + "/kmerClust_DONE"))

rule isONclustCon:
    input: OUTPUT_DIR + "/isONclustCon.fna"
    output: temp(touch(OUTPUT_DIR + "/isONclustCon_DONE"))

rule quant:
    input: OUTPUT_DIR + "/count_matrix.tsv"
    output: temp(touch(OUTPUT_DIR + "/quant_DONE"))

rule phylo: 
    input:
        expand(OUTPUT_DIR + "/sepp/ftable_filtered.{ext}", ext=["tsv", "biom"]),
        OUTPUT_DIR + "/sepp/sepp.tre",
    output: temp(touch(OUTPUT_DIR + "/phylo_DONE"))
