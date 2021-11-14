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

# call variants with medaka

# combine sample-specific amplicons and mapping