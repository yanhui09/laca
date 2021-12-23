checkpoint NGSpeciesID:
    input: rules.split_by_cluster.output
    output: directory(OUTPUT_DIR + "/NGSpeciesID/{barcode}/{c}")
    conda: "../envs/NGSpeciesID.yaml"
    params:
        racon_iter = config["NGSpeciesID"]["racon_iter"],
    log: OUTPUT_DIR + "/logs/NGSpeciesID/{barcode}/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/NGSpeciesID/{barcode}/{c}.txt"
    threads: config["threads"]["normal"]
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
    message: "Generate NGSpeciesID consensus (ID={wildcards.id}) in {wildcards.c} with medaka [{wildcards.barcode}]"
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
def get_isONclustCon(wildcards, pooling = True):
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_demultiplexed(wildcards)

    fnas = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            ids = glob_wildcards(checkpoints.NGSpeciesID.get(barcode=i, c=j).output[0] + "/reads_to_consensus_{id}.fastq").id
            for k in ids:
                fnas.append(OUTPUT_DIR + "/medaka/{barcode}/{c}/id_{id}/consensus.fasta".format(barcode=i, c=j, id=k))
    return fnas

rule collect_isONclustCon:
    input: lambda wc: get_isONclustCon(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/isONclustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i, id_i = i.split("/")[-4:-1]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "\n"
                        out.write(line)
