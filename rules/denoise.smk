# align merged assemblies with raw reads
# reused in racon iterations
rule minimap:
    input: 
      ref = OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/{assembly}.fna",
      fastq = OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    output: OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/{assembly}.paf",
    message: "{wildcards.c} for polish: alignments against {wildcards.assembly} assembly [{wildcards.barcode}]"
    params:
        x=config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/polish/minimap/{c}/{assembly}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/{barcode}/polish/minimap/{c}/{assembly}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/raw"
        return(prefix + ".paf", prefix + ".fna")
    else:
        prefix = OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/racon_{iter}".format(barcode=wildcards.barcode, c=wildcards.c, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
        get_racon_input,
    output: OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/racon_{iter}.fna"
    message: "Polish {wildcards.c} draft with racon, round={wildcards.iter} [{wildcards.barcode}]"
    params:
        m=config["racon"]["m"],
        x=config["racon"]["x"],
        g=config["racon"]["g"],
        w=config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR +"/logs/umap/{barcode}/polish/racon/{c}/round{iter}.log"
    benchmark: OUTPUT_DIR +"/benchmarks/umap/{barcode}/polish/racon/{c}/round{iter}.txt"
    threads: config["threads"]["large"]
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

checkpoint medaka_consensus:
    input:
        fna = expand(OUTPUT_DIR + "/umap/{{barcode}}/polish/draft/{{c}}/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    output: 
        fna = OUTPUT_DIR + "/umap/{barcode}/denoised_seqs/{c}.fna",
    message: "Generate consensus for {wildcards.c} with medaka [{wildcards.barcode}]"
    params:
        m=config["medaka"]["m"],
        _dir=OUTPUT_DIR + "/umap/{barcode}/polish/medaka/{c}",
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/polish/medaka/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/{barcode}/polish/medaka/{c}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fna} -o {params._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        cp {params._dir}/consensus.fasta {output.fna}
        """

# get {barcode} {c} from chekckpoint
def get_denoised_seqs(wildcards):
    barcodes = glob_wildcards(checkpoints.demultiplex.get(**wildcards).output[0]
     + "/{barcode, BRK[0-9][0-9]}/{runid}.fastq").barcode
    fnas = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.split_by_cluster.get(barcode=i).output[0] + "/{c}.fastq").c
        for j in cs:
            fnas.append(OUTPUT_DIR + "/umap/{barcode}/denoised_seqs/{c}.fna".format(barcode=i, c=j))
    return fnas

rule collect_denoised_seqs:
    input: get_denoised_seqs
    output: OUTPUT_DIR + "/umap/denoised_seqs.fna",
    message: "Collect denoised sequences"
    log: OUTPUT_DIR + "/logs/umap/denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/denoised_seqs.txt"
    shell:
        "cat {input} > {output} 2> {log}"

# to do: call variants with medaka

# dereplicate denoised sequences with mmseqs
rule dereplicate_denoised_seqs:
    input: OUTPUT_DIR + "/umap/denoised_seqs.fna",
    output: OUTPUT_DIR + "/rep_seqs.fna",
    message: "Dereplicate denoised sequences"
    log: OUTPUT_DIR + "/logs/dereplicate_denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/dereplicate_denoised_seqs.txt"
    params:
        tmp=OUTPUT_DIR + "/tmp",
    threads: config["threads"]["large"]
    conda: "../envs/mmseqs2.yaml"
    shell:
        "mmseqs easy-cluster {input} {output} {params.tmp} --threads {threads} --min-seq-id 1 -c 1 > {log} 2>&1"

# create abundance matrix with minimap2
