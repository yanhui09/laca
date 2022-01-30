rule isONclust:
    input: get_fq4Con(config["kmerbin"])
    output:
        _dir = directory(OUTPUT_DIR + "/isONclustCon/{barcode}/{c}"),
        tsv = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/final_clusters.tsv"
    conda: "../envs/isONcorCon.yaml"
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/isONclust.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/isONclust.txt"
    threads: config["threads"]["large"]
    shell:
        "isONclust --ont --fastq {input} "
        "--outfolder {output._dir} --t {threads} > {log} 2>&1"

def get_isONclust(wildcards, pooling = True, kmerbin = True):
    check_val("pooling", pooling, bool)
    check_val("kmerbin", kmerbin, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_demultiplexed(wildcards)

    bin2clusters = []
    for i in barcodes:
        if kmerbin == True:
            cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        else:
            cs = ["all"]
        for c in cs:
            bin2clusters.append(OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/final_clusters.tsv".format(barcode=i, c=c))
    return bin2clusters
 
checkpoint validate_clusters:
    input: lambda wc: get_isONclust(wc, pooling = config["pooling"], kmerbin = config["kmerbin"])
    output: directory(OUTPUT_DIR + "/isONclustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            barcode, c = [ i.split('/')[index] for index in [-3, -2] ]
            df_i = pd.read_csv(i, sep = '\t', header = None)
            df_i.columns = ['cluster', 'seqid']
            for clust_id, df_clust in df_i.groupby('cluster'):
                if len(df_clust) >= params.min_size:
                    df_clust['seqid'].to_csv(output[0] + "/{barcode}_{c}_{clust_id}.csv".format(barcode=barcode, c=c, clust_id=clust_id),
                     header = False, index = False)

rule get_fqs_split:
    input:
        bin2clust = OUTPUT_DIR + "/isONclustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        binned = get_fq4Con(config["kmerbin"]),
    output: OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/split/{clust_id}.fastq",
    conda: '../envs/seqkit.yaml'
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/{clust_id}/get_fqs_split.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/get_fqs_split.txt"
    shell: "seqkit grep -f {input.bin2clust} {input.binned} -o {output} --quiet 2> {log}"

rule spoa_consensus:
    input: rules.get_fqs_split.output
    output: OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/polish/{clust_id}/draft/raw.fna"
    conda: '../envs/spoa.yaml'
    params:
        l = config["spoa"]["l"],
        r = config["spoa"]["r"],
        g = config["spoa"]["g"],
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/{clust_id}/spoa_consensus.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/spoa_consensus.txt"
    shell: "spoa {input} -l {params.l} -r {params.r} -g {params.g} > {output} 2> {log}"

# polish with racon and medaka
use rule minimap2polish as minimap2polish2 with:
    input: 
        ref = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.fna",
        fastq = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/split/{clust_id}.fastq",
    output: 
        OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.paf",
    log: 
        OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/isOnclustCon/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.txt"

use rule racon as racon2 with:
    input:
        OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/split/{clust_id}.fastq",
        lambda wc: get_racon_input(wc, "isONclustCon"),
    output: 
        OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}.fna"
    log: 
        OUTPUT_DIR +"/logs/isONclustCon/{barcode}/{c}/{clust_id}/racon/round{iter}.log"
    benchmark: 
        OUTPUT_DIR +"/benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/racon/round{iter}.txt"

use rule medaka_consensus as medaka_consensus2 with:
    input:
        fna = expand(OUTPUT_DIR + "/isONclustCon/{{barcode}}/{{c}}/polish/{{clust_id}}/draft/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/split/{clust_id}.fastq",
    output: 
        fasta = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/polish/{clust_id}/medaka/consensus.fasta",
        _dir = directory(OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/polish/{clust_id}/medaka"),
    log: 
        OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/{clust_id}/medaka.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/medaka.txt"

# get {barcode} {c} {id} from chekckpoint
def get_isONclustCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.validate_clusters.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    fnas = []
    for i in bc_kb_cis:
        bc, kb, ci = i.split("_")
        fnas.append(OUTPUT_DIR + "/isONclustCon/{bc}/{kb}/polish/{ci}/medaka/consensus.fasta".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_isONclustCon:
    input: lambda wc: get_isONclustCon(wc),
    output: OUTPUT_DIR + "/isONclustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                bc_i, kb_i, ci_i = [ i.split("/")[index] for index in [-6, -5, -3] ]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + bc_i + "_" + kb_i + "_" + ci_i + "\n"
                        out.write(line)
