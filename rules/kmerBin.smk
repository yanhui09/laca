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

def get_bin(wildcards, pooling = True):
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_qced(wildcards)
    return expand(OUTPUT_DIR + "/kmerBin/{barcode}/hdbscan.tsv", barcode=barcodes)

# split reads by kmerbin
checkpoint cls_kmerbin:
    input: lambda wildcards: get_bin(wildcards, pooling = config["pooling"])
    output: directory(OUTPUT_DIR + "/kmerBin/clusters"),
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            bc = i.split("/")[-2]
            df = pd.read_csv(i, sep="\t")
            df = df[["read", "bin_id"]]
            kbs = df.bin_id.max()
            for kb in range(0, kbs+1):
                bc_kb = "/{bc}_c{kb}.csv".format(bc=bc, kb=kb)
                df.loc[df.bin_id == kb, "read"].to_csv(output[0] + bc_kb, header=False, index=False)
        
rule split_bin:
    input: 
        cluster = OUTPUT_DIR + "/kmerBin/clusters/{barcode}_{c}.csv",
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
    
    if kmerbin == True:
        fqs = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            fqs.append(OUTPUT_DIR + "/kmerBin/{bc}/split/{kb}.fastq".format(bc=bc, kb=kb))
    else:
        if pooling == True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        fqs = expand(OUTPUT_DIR + "/kmerBin/{bc}/all.fastq", bc=bcs)
    return fqs

def get_fq4Con(kmerbin = True):
    check_val("kmerbin", kmerbin, bool)
    if kmerbin == True:
        out = rules.split_bin.output
    else:
        out = rules.skip_bin.output
    return out