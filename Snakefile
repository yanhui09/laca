#--------------
configfile: "config.yaml"
#--------------
INPUT_DIR = config["basecalled_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
DATABASE_DIR = config["database_dir"].rstrip("/")

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

def check_val(var, val, class_type):
    if not isinstance(val, class_type):
        warns = ('\n\t' + str(var) + ' only accepts ' + str(class_type) + ' values.' +
         '\n\t' + str(val) + ' is used in config.yaml file.')
        raise ValueError(warns)

rule all:
    input:
        OUTPUT_DIR + "/demultiplex_DONE",
        OUTPUT_DIR + "/qc_DONE",
        OUTPUT_DIR + "/kmerClust_DONE",
        OUTPUT_DIR + "/isONclustCon_DONE",
        OUTPUT_DIR + "/isONcorCon_DONE",
        OUTPUT_DIR + "/quant_DONE",
        OUTPUT_DIR + "/taxonomy_DONE",

include: "rules/demultiplex.smk"
#include: "rules/umi.smk"
include: "rules/qc.smk"
include: "rules/kmerClust.smk"
include: "rules/isONclustCon.smk"
include: "rules/isONcorCon.smk"
include: "rules/quant.smk"
include: "rules/taxonomy.smk"
include: "rules/tree.smk"

rule demultiplex:
    input: lambda wc: expand(OUTPUT_DIR + "/raw/{barcode}.fastq", barcode=get_demultiplexed(wc))
    output: temp(touch(OUTPUT_DIR + "/.demultiplex_DONE"))

rule qc:
    input: lambda wc: get_filt(wc, pooling = config["pooling"]),
    output: temp(touch(OUTPUT_DIR + "/.qc_DONE"))

rule kmerClust:
    input: lambda wc: get_kmerClust(wc, pooling = config["pooling"]),
    output: temp(touch(OUTPUT_DIR + "/.kmerClust_DONE"))

rule isONclustCon:
    input: OUTPUT_DIR + "/isONclustCon.fna"
    output: temp(touch(OUTPUT_DIR + "/.isONclustCon_DONE"))

rule isONcorCon:
    input: OUTPUT_DIR + "/isONcorCon.fna"
    output: temp(touch(OUTPUT_DIR + "/.isONcorCon_DONE"))

rule quant:
    input: OUTPUT_DIR + "/count_matrix.tsv"
    output: temp(touch(OUTPUT_DIR + "/.quant_DONE"))

rule taxa: 
    input:
        [OUTPUT_DIR + "/" + str(x) + "/taxonomy.tsv" for x in config["classifier"]],
        OUTPUT_DIR + "/taxonomy.tsv",
    output: temp(touch(OUTPUT_DIR + "/.taxa_DONE"))

rule tree: 
    input:
        [OUTPUT_DIR + "/tree/" + str(x) + "/rooted_tree.nwk" for x in config["phylogen"]],
        OUTPUT_DIR + "/tree.nwk",
    output: temp(touch(OUTPUT_DIR + "/.tree_DONE"))
