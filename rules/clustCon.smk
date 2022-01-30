# draw draft with max average score from pairwise alignments
rule minimap2clust:
    input: get_fq4Con(config["kmerbin"])
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
    output: OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust.csv",
    conda: '../envs/clustCon.yaml'
    params:
        prefix = OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust",
        min_score_frac = config["clustCon"]["min_score_frac"],
        min_reads = config["min_cluster_size"],
        max_recurs = config["clustCon"]["max_recursion"],
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/bin2clust/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/bin2clust/{c}.txt"
    shell:
        "python scripts/binClust.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input} > {log} 2>& 1"

def get_clust(wildcards, pooling = True, kmerbin = True):
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
            bin2clusters.append(OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust.csv".format(barcode=i, c=c))
    return bin2clusters
 
checkpoint clusters_cleanup:
    input: lambda wc: get_clust(wc, pooling = config["pooling"], kmerbin = config["kmerbin"])
    output: directory(OUTPUT_DIR + "/clustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import shutil
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            num_lines = sum(1 for line in open(i))
            if num_lines > params.min_size:
                barcode, c = [ i.split("/")[index] for index in [-4, -2] ]
                #shutil.copy(i, output[0] + "/{barcode}_{c}.csv".format(barcode=barcode, c=c))
                df = pd.read_csv(i)
                for clust_id, df_clust in df.groupby('cluster'):
                    if len(df_clust) >= params.min_size:
                        bc_kb_ci = "/{barcode}_{c}_id{clust_id}".format(barcode=barcode, c=c, clust_id=clust_id)
                        df_clust['read_id'].to_csv(output[0] + bc_kb_ci + ".csv", header = False, index = False)
                        ref_idx  = df_clust['clust_read_score'].idxmax()
                        ref_read = df_clust.loc[ref_idx, ['read_id']]
                        ref_read.to_csv(output[0] + bc_kb_ci + ".ref", header=False, index=False)
   
rule get_fqs_split2:
    input:
        pool = OUTPUT_DIR + "/clustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        ref = OUTPUT_DIR + "/clustCon/clusters/{barcode}_{c}_{clust_id}.ref",
        binned = get_fq4Con(config["kmerbin"]),
    output:
        pool = OUTPUT_DIR + "/clustCon/{barcode}/{c}/split/{clust_id}.fastq",
        ref = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/raw.fna",
    conda: '../envs/seqkit.yaml'
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/{c}/{clust_id}/get_fqs_split2.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/{c}/{clust_id}/get_fqs_split2.txt"
    shell:
        """
        seqkit grep -f {input.pool} {input.binned} -o {output.pool} --quiet 2> {log}
        seqkit grep -f {input.ref} {input.binned} --quiet | seqkit fq2fa | \
        seqkit replace -p '^(.+)$' -r '{wildcards.barcode}_{wildcards.c}_{wildcards.clust_id}' -o {output.ref} 2> {log}
        """

# align assemblies with raw reads
# reused in racon iterations
rule minimap2polish:
    input: 
      ref = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.fna",
      fastq = OUTPUT_DIR + "/clustCon/{barcode}/{c}/split/{clust_id}.fastq",
    output: OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.paf",
    message: "Polish {wildcards.c} draft [id={wildcards.clust_id}]: alignments against {wildcards.assembly} assembly [{wildcards.barcode}]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.txt"
    threads: config["threads"]["normal"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards, cluster):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUTPUT_DIR + "/" + cluster + "/{barcode}/{c}/polish/{clust_id}/draft/raw"
        return(prefix + ".paf", prefix + ".fna")
    else:
        prefix = OUTPUT_DIR + "/" + cluster + "/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}".format(barcode=wildcards.barcode,
         c=wildcards.c, clust_id=wildcards.clust_id, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        OUTPUT_DIR + "/clustCon/{barcode}/{c}/split/{clust_id}.fastq",
        lambda wc: get_racon_input(wc, "clustCon"),
    output: OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}.fna"
    message: "Polish {wildcards.c} draft [id={wildcards.clust_id}] with racon, round={wildcards.iter} [{wildcards.barcode}]"
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR +"/logs/clustCon/{barcode}/{c}/{clust_id}/racon/round{iter}.log"
    benchmark: OUTPUT_DIR +"/benchmarks/clustCon/{barcode}/{c}/{clust_id}/racon/round{iter}.txt"
    threads: config["threads"]["normal"]
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

rule medaka_consensus:
    input:
        fna = expand(OUTPUT_DIR + "/clustCon/{{barcode}}/{{c}}/polish/{{clust_id}}/draft/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = OUTPUT_DIR + "/clustCon/{barcode}/{c}/split/{clust_id}.fastq",
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
        -d {input.fna} -o {output._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        """

# get {barcode} {c} from chekckpoint
def get_clustCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.clusters_cleanup.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    #b_cs = glob_wildcards(OUTPUT_DIR + "/clustCon/clusters/{b_c}.csv").b_c
    fnas = []
    for i in bc_kb_cis:
        bc, kb, ci = i.split("_")
        fnas.append(OUTPUT_DIR + "/clustCon/{bc}/{kb}/polish/{ci}/medaka/consensus.fasta".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_clustCon:
    input: lambda wc: get_clustCon(wc),
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
