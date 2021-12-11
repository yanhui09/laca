# linker and primer info
fprimers = config["fprimer"]
rprimers = config["rprimer"]

# reverse complementation
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
fprimersR = {k: revcomp(v) for k, v in fprimers.items()}
rprimersR = {k: revcomp(v) for k, v in rprimers.items()}

# pattern search for umi using cutadapt
# nanopore possibly sequences either strand
def seqs_join(primer1, primer2):
    joined = '-g ' + primer1 + '...' + revcomp(primer2)
    return joined
def linked_pattern(primers1, primers2):
    primers1_values = list(primers1.values())
    primers2_values = list(primers2.values())
    linked = [seqs_join(primer1, primer2) for primer1 in primers1_values for primer2 in primers2_values]
    return ' '.join(linked)

# pattern
f5_pattern1 = linked_pattern(fprimers, rprimers)
f5_pattern2 = linked_pattern(rprimers, fprimers)
f5_patterns = f5_pattern1 + ' ' + f5_pattern2
#---------

rule combine_fastq:
    input: get_demultiplexed
    output: temp(OUTPUT_DIR + "/raw_fq/pooled.fastq")
    shell:
        "cat {input} > {output}"

# trim primers 
rule trim_primers:
    input: OUTPUT_DIR + "/raw_fq/{barcode}.fastq"
    output: OUTPUT_DIR + "/umap/{barcode}/trimmed.fastq",
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5_patterns,
        max_err = config["max_err"],
        min_overlap = config["min_overlap"],
    log: OUTPUT_DIR + "/logs/umap/{barcode}/trim_primers.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/{barcode}/trim_primers.txt"
    threads: config["threads"]["large"]
    shell:
        """
        cutadapt \
            -j {threads} -e {params.max_err} -O {params.min_overlap} \
            --discard-untrimmed \
            {params.f} \
            -o {output} \
            {input} \
            > {log} 2>&1
        """

# quality filter
rule nanofilt_umap:
    input:  rules.trim_primers.output,
    output: OUTPUT_DIR + "/umap/{barcode}/filt.fastq"
    conda: "../envs/nanofilt.yaml"
    params:
        q = config["NanoFilt"]["q"],
        l = config["NanoFilt"]["l"],
        maxlength = config["NanoFilt"]["maxlength"],
    log: OUTPUT_DIR + "/logs/umap/{barcode}/nanofilt.log"
    shell: 
        "cat {input} | NanoFilt -q {params.q} -l {params.l} --maxlength {params.maxlength} 2> {log} 1> {output}"

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

rule checkend_umap_hdbscan:
    input: lambda wc: get_kmerClust(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/KmerClust_DONE",
    shell: "touch {output}"