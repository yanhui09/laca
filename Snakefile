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

# check list elements
def check_list_ele(var, val, eles):
    var = var.capitalize()
    if val:
        for i in val:
            if i not in eles:
                raise ValueError("\t\n{} parameters not recognized.\t\nPlease choose from {} in the config.yaml file.".format(var,eles))
    else:
        raise ValueError("\t\n{} parameters not specified.\t\nPlease choose from {} in the config.yaml file.".format(var,eles))

rule all:
    input:
        OUTPUT_DIR + "/.demultiplex_DONE",
        OUTPUT_DIR + "/.qc_DONE",
        OUTPUT_DIR + "/.kmerBin_DONE",
        OUTPUT_DIR + "/.quant_DONE",
        OUTPUT_DIR + "/.taxa_DONE",
        OUTPUT_DIR + "/.tree_DONE",

include: "rules/init.smk"
#include: "rules/nanosim.smk"
include: "rules/demultiplex.smk"
include: "rules/qc.smk"
include: "rules/kmerBin.smk"
include: "rules/clustCon.smk"
include: "rules/umiCon.smk"
include: "rules/quant.smk"
include: "rules/requant.smk"
include: "rules/taxonomy.smk"
include: "rules/tree.smk"

rule demultiplex:
    input: lambda wc: expand(OUTPUT_DIR + "/qc/{barcode}.fastq", barcode=get_demultiplexed(wc))
    output: temp(touch(OUTPUT_DIR + "/.demultiplex_DONE"))

rule qc:
    input: lambda wc: get_filt(wc, pooling = config["pooling"]),
    output: temp(touch(OUTPUT_DIR + "/.qc_DONE"))

rule kmerBin:
    input: lambda wc: get_kmerBin(wc, pooling = config["pooling"], kmerbin = config["kmerbin"]),
    output: temp(touch(OUTPUT_DIR + "/.kmerBin_DONE"))

rule kmerCon:
    input:
        OUTPUT_DIR + "/.kmerBin_DONE",
        OUTPUT_DIR + "/kmerCon.fna"
    output: temp(touch(OUTPUT_DIR + "/.kmerCon_DONE"))

rule clustCon:
    input:
        OUTPUT_DIR + "/.kmerBin_DONE",
        OUTPUT_DIR + "/clustCon.fna"
    output: temp(touch(OUTPUT_DIR + "/.clustCon_DONE"))

rule isONclustCon:
    input:
        OUTPUT_DIR + "/.kmerBin_DONE",
        OUTPUT_DIR + "/isONclustCon.fna"
    output: temp(touch(OUTPUT_DIR + "/.isONclustCon_DONE"))

rule isONcorCon:
    input: 
        OUTPUT_DIR + "/.kmerBin_DONE",
        OUTPUT_DIR + "/isONcorCon.fna"
    output: temp(touch(OUTPUT_DIR + "/.isONcorCon_DONE"))

rule umiCon:
    input: 
        OUTPUT_DIR + "/umiCon.fna"
    output: temp(touch(OUTPUT_DIR + "/.umiCon_DONE"))

rule quant:
    input: chimeraF(config["chimeraF"])
    output: temp(touch(OUTPUT_DIR + "/.quant_DONE"))

rule taxa: 
    input:
        [OUTPUT_DIR + "/taxonomy/" + str(x) + "/taxonomy.tsv" for x in config["classifier"]],
        OUTPUT_DIR + "/taxonomy.tsv",
    output: temp(touch(OUTPUT_DIR + "/.taxa_DONE"))

rule tree: 
    input:
        [OUTPUT_DIR + "/tree/" + str(x) + "/tree.nwk" for x in config["phylogen"]],
        OUTPUT_DIR + "/tree.nwk",
    output: temp(touch(OUTPUT_DIR + "/.tree_DONE"))

rule requant:
    input:
        OUTPUT_DIR + "/rep_seqs_requant.fasta",
        OUTPUT_DIR + "/count_matrix_requant.tsv",
    output: temp(touch(OUTPUT_DIR + "/.requant_DONE"))
