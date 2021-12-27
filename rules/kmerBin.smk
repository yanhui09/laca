# kmer calculation
rule kmer_freqs:
    input: OUTPUT_DIR + "/qc/qfilt/{barcode}.fastq"
    output: OUTPUT_DIR + "/kmerBin/{barcode}/kmer_freqs.txt"
    conda: "../envs/kmerBin.yaml"
    params: 
        kmer_size = config["kmer_size"],
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/kmer_freqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/kmer_freqs.txt"
    threads: config["threads"]["large"]
    shell:
        "python scripts/kmerFreqs.py"
        " -k {params.kmer_size}"
        " -r {input} -t {threads}"
        " 2> {log} > {output}"

# kmer binning
rule umap:
    input: rules.kmer_freqs.output
    output: 
        cluster=OUTPUT_DIR + "/kmerBin/{barcode}/hdbscan.tsv",
	    plot=OUTPUT_DIR + "/kmerBin/{barcode}/hdbscan.png",
    conda: "../envs/kmerBin.yaml"
    params:
        n_neighbors = config["umap"]["n_neighbors"],
        min_dist = config["umap"]["min_dist"],
        n_components = config["umap"]["n_components"],
	    min_cluster_size = config["hdbscan"]["min_cluster_size"],
        min_samples = config["hdbscan"]["min_samples"],
	    epsilon = config["hdbscan"]["epsilon"],
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/umap.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/umap.txt"
    threads: config["threads"]["large"]
    shell:
       "NUMBA_NUM_THREADS={threads} python scripts/kmerBin.py -k {input}"
       " -n {params.n_neighbors} -d {params.min_dist} -t {params.n_components}"
       " -s {params.min_cluster_size} -m {params.min_samples} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

# split reads by cluster
checkpoint cluster_info:
    input: rules.umap.output.cluster,
    output: directory(OUTPUT_DIR + "/kmerBin/{barcode}/clusters"),
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/clusters.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/clusters.txt"
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
        clusters = OUTPUT_DIR + "/kmerBin/{barcode}/clusters/{c}.txt",
        fastq = OUTPUT_DIR + "/raw/qfilt/{barcode}.fastq",
    output: OUTPUT_DIR + "/kmerBin/{barcode}/clusters/{c}.fastq",
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/clusters/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/clusters/{c}.txt"
    shell:
        "seqkit grep {input.fastq} -f {input.clusters} -o {output} 2> {log}"

# get {barcode} {c} from chekckpoint
def get_kmerBin(wildcards, pooling = True):
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_demultiplexed(wildcards)

    fqs = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            fqs.append(OUTPUT_DIR + "/kmerBin/{barcode}/clusters/{c}.fastq".format(barcode=i, c=j))
    return fqs
