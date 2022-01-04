checkpoint NGSpeciesID:
    input: get_fq4Con(config["kmerbin"])
    output: directory(OUTPUT_DIR + "/isONclustCon/{barcode}/{c}")
    conda: "../envs/NGSpeciesID.yaml"
    params:
        racon_iter = config["NGSpeciesID"]["racon_iter"],
        abundance_ratio = config["NGSpeciesID"]["abundance_ratio"],
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/NGSpeciesID.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/NGSpeciesID.txt"
    threads: config["threads"]["normal"]
    shell:
        "NGSpeciesID --ont --consensus"
        " --racon --racon_iter {params.racon_iter}"
        " --abundance_ratio {params.abundance_ratio}"
        " --fastq {input} --outfolder {output} --t {threads} > {log} 2>&1"

rule medaka:
    input:
        fasta = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/racon_cl_id_{id}/consensus.fasta",
        fastq = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/reads_to_consensus_{id}.fastq",
    output: 
        fasta = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/id_{id}/medaka/consensus.fasta",
        _dir = directory(OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/id_{id}/medaka"),
    message: "Generate NGSpeciesID consensus (ID={wildcards.id}) in {wildcards.c} with medaka [{wildcards.barcode}]"
    params:
        m = config["medaka"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/medaka/id_{id}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/medaka/id_{id}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fasta} -o {output._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        """

# get {barcode} {c} {id} from chekckpoint
def get_isONclustCon(wildcards, pooling = True, kmerbin = True):
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_demultiplexed(wildcards)

    fnas = []
    for i in barcodes:
        if kmerbin == True:
            cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        else:
            cs = ["all"]
        for j in cs:
            ids = glob_wildcards(checkpoints.NGSpeciesID.get(barcode=i, c=j).output[0] + "/reads_to_consensus_{id}.fastq").id
            for k in ids:
                fnas.append(OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/id_{id}/medaka/consensus.fasta".format(barcode=i, c=j, id=k))
    return fnas

rule collect_isONclustCon:
    input: lambda wc: get_isONclustCon(wc, pooling = config["pooling"], kmerbin = config["kmerbin"]),
    output: OUTPUT_DIR + "/isONclustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i, id_i = i.split("/")[-5:-2]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "\n"
                        out.write(line)
