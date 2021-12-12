# kmer calculation
rule kmer_freqs:
    input: rules.nanofilt_umap.output
    output: OUTPUT_DIR + "/umap/{barcode}/kmer_freqs.txt"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/kmer_freqs.log"
    threads: config["threads"]["normal"]
    conda: "../envs/kmer_freqs.yaml"
    params: 
        scripts = "scripts",
        kmer_size = config["kmer_size"],
    shell:
        "python {params.scripts}/kmer_freqs.py"
        " -k {params.kmer_size}"
        " -r {input} -t {threads}"
        " 2> {log} > {output}"

# umap cluster
rule umap:
    input: rules.kmer_freqs.output
    output: 
        cluster=OUTPUT_DIR + "/umap/{barcode}/hdbscan.tsv",
	    plot=OUTPUT_DIR + "/umap/{barcode}/hdbscan.png",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/umap.log"
    conda: "../envs/umap_cluster.yaml"
    params:
        n_neighbors = config["umap"]["n_neighbors"],
        min_dist = config["umap"]["min_dist"],
        n_components = config["umap"]["n_components"],
	    min_cluster_size = config["hdbscan"]["min_cluster_size"],
        min_samples = config["hdbscan"]["min_samples"],
	    epsilon = config["hdbscan"]["epsilon"],
    threads: config["threads"]["large"]
    shell:
       "NUMBA_NUM_THREADS={threads} python scripts/umap_cluster.py -k {input}"
       " -n {params.n_neighbors} -d {params.min_dist} -t {params.n_components}"
       " -s {params.min_cluster_size} -m {params.min_samples} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

# split reads by cluster
checkpoint cluster_info:
    input: rules.umap.output.cluster,
    output: directory(OUTPUT_DIR + "/umap/{barcode}/clusters"),
    log: OUTPUT_DIR + "/logs/umap/{barcode}/clusters.log"
    params:
        scripts = "scripts"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        df = df[["read", "bin_id"]]
        clusters = df.bin_id.max()
        os.makedirs(output[0])
        for cluster in range(0, clusters+1):
            df.loc[df.bin_id == cluster, "read"].to_csv(output[0] + "/c" + str(cluster) + ".txt", sep="\t", header=False, index=False)
        
rule split_by_cluster:
    input: 
        clusters = OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.txt",
        fastq = rules.nanofilt_umap.output,
    output: OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/clusters/{c}.log"
    conda: "../envs/seqkit.yaml"
    shell:
        "seqkit grep {input.fastq} -f {input.clusters} -o {output} 2> {log}"

# get {barcode} {c} from chekckpoint
def get_kmerClust(wildcards, pooling = True):
    if pooling:
        barcodes = ["pooled"]
    elif pooling is False:
        barcodes = glob_wildcards(checkpoints.guppy_demultiplex.get(**wildcards).output[0]
        + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    else:
        raise ValueError('Pooling only allows bool type [True/False].\n{} is used in the config file'.format(x))

    fqs = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            fqs.append(OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq".format(barcode=i, c=j))
    return fqs
