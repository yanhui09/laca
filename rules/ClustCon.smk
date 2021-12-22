# draw draft with max average score from pairwise alignments
rule minimap2clust:
    input: rules.split_by_cluster.output
    output: OUTPUT_DIR + "/ClustCon/{barcode}/avr_aln/{c}/minimap2clust.paf"
    conda: '../envs/polish.yml'
    log: OUTPUT_DIR + "/ClustCon/{barcode}/minimap2clust/{c}.log"
    benchmark: OUTPUT_DIR + "/ClustCon/{barcode}/minimap2clust/{c}.text"
    shell:
        "minimap2 -t {threads} -x ava-ont --no-long-join -r100"
        " {input} {input} > {output} 2> {log}"

rule bin2clustering:
    input: OUTPUT_DIR + "/ClustCon/{barcode}/avr_aln/{c}/minimap2clust.paf"
    output:
        heatmap = OUTPUT_DIR + "/ClustCon/{barcode}/avr_aln/{c}/bin2clust.heatmap.png",
        info = OUTPUT_DIR + "/ClustCon/{barcode}/avr_aln/{c}/bin2clust.info.csv",
    conda: '../envs/ClustCon.yaml'
    params:
        prefix = OUTPUT_DIR + "/ClustCon/{barcode}/avr_aln/{c}/bin2clust",
        min_score_frac = config["ClustCon"]["min_score_frac"],
        min_reads = config["ClustCon"]["min_reads"],
        max_recurs = config["ClustCon"]["max_recursion"],
    shell:
        "python scripts/cluster_ava_alignments.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input}"


# align merged assemblies with raw reads
# reused in racon iterations
rule minimap2polish:
    input: 
      ref = OUTPUT_DIR + "/ClustCon/{barcode}/polish/draft/{c}/{assembly}.fna",
      fastq = OUTPUT_DIR + "/ClustCon/{barcode}/clusters/{c}.fastq",
    output: OUTPUT_DIR + "/ClustCon/{barcode}/polish/draft/{c}/{assembly}.paf",
    message: "{wildcards.c} for polish: alignments against {wildcards.assembly} assembly [{wildcards.barcode}]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/ClustCon/{barcode}/polish/minimap/{c}/{assembly}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/ClustCon/{barcode}/polish/minimap/{c}/{assembly}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUTPUT_DIR + "/ClustCon/{barcode}/polish/draft/{c}/raw"
        return(prefix + ".paf", prefix + ".fna")
    else:
        prefix = OUTPUT_DIR + "/ClustCon/{barcode}/polish/draft/{c}/racon_{iter}".format(barcode=wildcards.barcode, c=wildcards.c, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        OUTPUT_DIR + "/ClustCon/{barcode}/clusters/{c}.fastq",
        get_racon_input,
    output: OUTPUT_DIR + "/ClustCon/{barcode}/polish/draft/{c}/racon_{iter}.fna"
    message: "Polish {wildcards.c} draft with racon, round={wildcards.iter} [{wildcards.barcode}]"
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR +"/logs/ClustCon/{barcode}/polish/racon/{c}/round{iter}.log"
    benchmark: OUTPUT_DIR +"/benchmarks/ClustCon/{barcode}/polish/racon/{c}/round{iter}.txt"
    threads: 1
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

checkpoint medaka_consensus:
    input:
        fna = expand(OUTPUT_DIR + "/ClustCon/{{barcode}}/polish/draft/{{c}}/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = OUTPUT_DIR + "/ClustCon/{barcode}/clusters/{c}.fastq",
    output: 
        fna = OUTPUT_DIR + "/ClustCon/{barcode}/denoised_seqs/{c}.fna",
    message: "Generate consensus for {wildcards.c} with medaka [{wildcards.barcode}]"
    params:
        m = config["medaka"]["m"],
        _dir = OUTPUT_DIR + "/ClustCon/{barcode}/polish/medaka/{c}",
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/ClustCon/{barcode}/polish/medaka/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/ClustCon/{barcode}/polish/medaka/{c}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fna} -o {params._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        cp {params._dir}/consensus.fasta {output.fna}
        """

# get {barcode} {c} from chekckpoint
def get_ClustCon(wildcards, pooling = True):
    if pooling:
        barcodes = ["pooled"]
    elif pooling is False:
        barcodes = glob_wildcards(checkpoints.demultiplex.get(**wildcards).output[0]
        + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    else:
        raise ValueError('Pooling only allows bool type [True/False].\n{} is used in the config file'.format(x))

    fnas = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            fnas.append(OUTPUT_DIR + "/ClustCon/{barcode}/denoised_seqs/{c}.fna".format(barcode=i, c=j))
    return fnas

rule collect_ClustCon:
    input: lambda wc: get_ClustCon(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/ClustCon.fna",
    message: "Collect consensus sequences from kmer clusters."
    log: OUTPUT_DIR + "/logs/ClustCon/collect_ClustCon.log"
    benchmark: OUTPUT_DIR + "/benchmarks/ClustCon/collect_ClustCon.txt"
    shell:
        "cat {input} > {output} 2> {log}"