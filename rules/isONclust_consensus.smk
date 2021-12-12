checkpoint NGSpeciesID:
    input: rules.split_by_cluster.output
    output: directory(OUTPUT_DIR + "/NGSpeciesID/{barcode}/{c}")
    threads: config["threads"]["normal"]
    conda: "../envs/NGSpeciesID.yaml"
    params:
        racon_iter = config["NGSpeciesID"]["racon_iter"],
    log: OUTPUT_DIR + "/logs/NGSpeciesID/{barcode}/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/NGSpeciesID/{barcode}/{c}.txt"
    shell:
        "NGSpeciesID --ont --consensus"
        " --racon --racon_iter {params.racon_iter}"
        " --fastq {input} --outfolder {output} --t {threads} > {log} 2>&1"

rule medaka:
    input:
        fasta = OUTPUT_DIR + "/NGSpeciesID/{barcode}/{c}/racon_cl_id_{id}/consensus.fasta",
        fastq = OUTPUT_DIR + "/NGSpeciesID/{barcode}/{c}/reads_to_consensus_{id}.fastq",
    output: 
        fasta = OUTPUT_DIR + "/medaka/{barcode}/{c}/id_{id}/consensus.fasta",
        _dir = directory(OUTPUT_DIR + "/medaka/{barcode}/{c}/id_{id}"),
    message: "Polish NGSpeciesID consensus (ID={wildcards.id}) in {wildcards.c} with medaka [{wildcards.barcode}]"
    params:
        m = config["medaka"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/medaka/{barcode}/{c}/id_{id}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/medaka/{barcode}/{c}/id_{id}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fasta} -o {output._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        """

# get {barcode} {c} {id} from chekckpoint
def get_isONclust_consensus(wildcards, pooling = True):
    if pooling:
        barcodes = ["pooled"]
    elif pooling is False:
        barcodes = glob_wildcards(checkpoints.guppy_demultiplex.get(**wildcards).output[0]
        + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
        barcodes = list(set(barcodes))
    else:
        raise ValueError('Pooling only allows bool type [True/False].\n{} is used in the config file'.format(x))

    fnas = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            ids = glob_wildcards(checkpoints.NGSpeciesID.get(barcode=i, c=j).output[0] + "/reads_to_consensus_{id}.fastq").id
            for k in ids:
                fnas.append(OUTPUT_DIR + "/medaka/{barcode}/{c}/id_{id}/consensus.fasta".format(barcode=i, c=j, id=k))
    return fnas

rule collect_isONclust_consensus:
    input: lambda wc: get_isONclust_consensus(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/isONclustConsensus.fna"
    shell: "cat {input} > {output}"
