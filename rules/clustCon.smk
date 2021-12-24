# draw draft with max average score from pairwise alignments
rule minimap2clust:
    input: rules.split_by_cluster.output
    output: OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/minimap2clust.paf"
    conda: '../envs/polish.yaml'
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/minimap2clust/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/minimap2clust/{c}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x ava-ont --no-long-join -r100"
        " {input} {input} > {output} 2> {log}"

rule bin2clustering:
    input: OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/minimap2clust.paf"
    output:
        heatmap = OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust.heatmap.png",
        info = OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust.info.csv",
    conda: '../envs/clustCon.yaml'
    params:
        prefix = OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust",
        min_score_frac = config["clustCon"]["min_score_frac"],
        min_reads = config["clustCon"]["min_reads"],
        max_recurs = config["clustCon"]["max_recursion"],
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/bin2clust/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/bin2clust/{c}.txt"
    shell:
        "python scripts/cluster_ava_alignments.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input} > {log} 2>& 1"

checkpoint get_bin2clust:
    input: rules.bin2clustering.output.info,
    output: directory(OUTPUT_DIR + "/clustCon/{barcode}/{c}/clusters")
    run:
        import pandas as pd
        df = pd.read_csv(input[0])
        for clust_id, df_clust in df.groupby('cluster'):
            clust_id = f"id_{clust_id}"
            clust_dir = output[0] + str("/{clust_id}").format(clust_id=clust_id)
            os.makedirs(clust_dir, exist_ok=True)
            read_pool = clust_dir + "/pool.csv"
            df_clust['read_id'].to_csv(read_pool, index=False, header=False)
            read_ref = clust_dir + "/ref.csv"
            ref_idx  = df_clust['clust_read_score'].idxmax()
            ref_read = df_clust.loc[ref_idx, ['read_id']]
            ref_read.to_csv(read_ref, index=False, header=False)
    
rule get_clust_reads:
    input:
        pool = OUTPUT_DIR + "/clustCon/{barcode}/{c}/clusters/{clust_id}/pool.csv",
        ref = OUTPUT_DIR + "/clustCon/{barcode}/{c}/clusters/{clust_id}/ref.csv",
        binned = rules.split_by_cluster.output,
    output:
        pool = OUTPUT_DIR + "/clustCon/{barcode}/{c}/clusters/{clust_id}/pool.fastq",
        ref = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/raw.fna",
    conda: '../envs/seqkit.yaml'
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/{c}/id_{clust_id}/get_clust_reads.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/{c}/id_{clust_id}/get_clust_reads.txt"
    shell:
        """
        seqkit grep -f {input.pool} {input.binned} > {output.pool} 2> {log}
        seqkit grep -f {input.ref} {input.binned} | seqkit fq2fa > {output.ref} 2> {log}
        """

# align merged assemblies with raw reads
# reused in racon iterations
rule minimap2polish:
    input: 
      ref = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.fna",
      fastq = OUTPUT_DIR + "/clustCon/{barcode}/{c}/clusters/{clust_id}/pool.fastq",
    output: OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.paf",
    message: "Polish {wildcards.c} [id={wildcards.clust_id}]: alignments against {wildcards.assembly} assembly [{wildcards.barcode}]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/{c}/id_{clust_id}/minimap2polish/{assembly}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/{c}/id_{clust_id}/minimap2polish/{assembly}.txt"
    threads: config["threads"]["normal"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/raw"
        return(prefix + ".paf", prefix + ".fna")
    else:
        prefix = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}".format(barcode=wildcards.barcode,
         c=wildcards.c, clust_id=wildcards.clust_id, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        OUTPUT_DIR + "/clustCon/{barcode}/{c}/clusters/{clust_id}/pool.fastq",
        get_racon_input,
    output: OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}.fna"
    message: "Polish {wildcards.c} draft [id={wildcards.clust_id}] with racon, round={wildcards.iter} [{wildcards.barcode}]"
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR +"/logs/clustCon/{barcode}/{c}/id_{clust_id}/racon/round{iter}.log"
    benchmark: OUTPUT_DIR +"/benchmarks/clustCon/{barcode}/{c}/id_{clust_id}/racon/round{iter}.txt"
    threads: 1
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

rule medaka_consensus:
    input:
        fna = expand(OUTPUT_DIR + "/clustCon/{{barcode}}/{{c}}/polish/{{clust_id}}/draft/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = OUTPUT_DIR + "/clustCon/{barcode}/{c}/clusters/{clust_id}/pool.fastq",
    output: 
        fasta = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/medaka/consensus.fasta",
        _dir = directory(OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/medaka"),
    message: "Generate Clust consensus [id={wildcards.clust_id}] in {wildcards.c} with medaka [{wildcards.barcode}]"
    params:
        m = config["medaka"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/{c}/id_{clust_id}/medaka.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/{c}/id_{clust_id}/medaka.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fasta} -o {output._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        """

# get {barcode} {c} from chekckpoint
def get_clustCon(wildcards, pooling = True):
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_demultiplexed(wildcards)

    fnas = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            ids = glob_wildcards(checkpoints.get_bin2clust.get(barcode=i, c=j).output[0] + "/id_{clust_id}/pool.csv").clust_id
            for k in ids:
                fnas.append(OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/id_{id}/medaka/consensus.fasta".format(barcode=i, c=j, id=k))
    return fnas

rule collect_clustCon:
    input: lambda wc: get_clustCon(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/clustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i, id_i = [ i.split("/")[index] for index in [-6, -5, -3] ]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "\n"
                        out.write(line)
