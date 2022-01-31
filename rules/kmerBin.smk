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
	    min_bin_size = config["hdbscan"]["min_bin_size"],
        min_samples = config["hdbscan"]["min_samples"],
	    epsilon = config["hdbscan"]["epsilon"],
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/umap.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/umap.txt"
    threads: config["threads"]["large"]
    shell:
       "NUMBA_NUM_THREADS={threads} python scripts/kmerBin.py -k {input}"
       " -n {params.n_neighbors} -d {params.min_dist} -t {params.n_components}"
       " -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

# split reads by cluster
checkpoint cls_kmerbin:
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
        
rule split_bin:
    input: 
        cluster = OUTPUT_DIR + "/kmerBin/{barcode}/clusters/{c}.txt",
        fqs = OUTPUT_DIR + "/qc/qfilt/{barcode}.fastq",
    output: OUTPUT_DIR + "/kmerBin/{barcode}/split/{c}.fastq",
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/split/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/split/{c}.txt"
    shell:
        "seqkit grep {input.fqs} -f {input.cluster} -o {output} 2> {log}"

rule skip_bin:
    input: OUTPUT_DIR + "/qc/qfilt/{barcode}.fastq"
    output: OUTPUT_DIR + "/kmerBin/{barcode}/all.fastq"
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/skip_bin.log"
    shell: "cp -p {input} {output} 2> {log}"

# get {barcode} {c} from chekckpoint
def get_kmerBin(wildcards, pooling = True, kmerbin = True):
    check_val("pooling", pooling, bool)
    check_val("kmerbin", kmerbin, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_demultiplexed(wildcards)

    fqs = []
    for i in barcodes:
        if kmerbin == True:
            cs = glob_wildcards(checkpoints.cls_kmerbin.get(barcode=i).output[0] + "/{c}.txt").c
            for j in cs:
                fqs.append(OUTPUT_DIR + "/kmerBin/{barcode}/split/{c}.fastq".format(barcode=i, c=j))
        else:
            fqs.append(OUTPUT_DIR + "/kmerBin/{barcode}/all.fastq".format(barcode=i))
    return fqs

def get_fq4Con(kmerbin = True):
    check_val("kmerbin", kmerbin, bool)
    if kmerbin == True:
        out = rules.split_bin.output
    else:
        out = rules.skip_bin.output
    return out