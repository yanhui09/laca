# kmer calculation
rule kmer_freqs:
    input: "qc/qfilt/{barcode}.fastq"
    output: "kmerBin/{barcode}/kmer_freqs.txt"
    conda: "../envs/kmerBin.yaml"
    params: 
        kmer_size = config["kmer_size"],
    log: "logs/kmerBin/{barcode}/kmer_freqs.log"
    benchmark: "benchmarks/kmerBin/{barcode}/kmer_freqs.txt"
    threads: config["threads"]["large"]
    shell:
        "python scripts/kmerFreqs.py"
        " -k {params.kmer_size}"
        " -r {input} -t {threads}"
        " 2> {log} > {output}"

def get_batch_size(batch_size, ram, kmer_file):
    # necessary in dry run
    if not os.path.exists(kmer_file):
        return -1 # pesduo
    # total input read size
    with open(kmer_file) as f:
        total = sum(1 for line in f) -1
    
    # zero and negtive values to use all reads
    try: 
        batch_size = int(batch_size)
        if batch_size <= 0:
            batch_size = total
        # auto estimation of batch size base on RAM
    except:
        if batch_size == "auto":
            # estimated constan m from 333 reads with max RAM of 241.26 MB
            m = 241.26 * 2 / 333 ** 2
            est_size = int((2 * 1024 * ram / m) ** 0.5)
            if est_size > total:
                batch_size = total
            else:
                batch_size = est_size
    return batch_size

# check whether to split the kmer freqs in batches
checkpoint shuffle_batch:
    input: rules.kmer_freqs.output
    output: directory("kmerBin/{barcode}/batch_check")
    conda: "../envs/coreutils.yaml"
    params:
        batch_size = lambda wc, input: get_batch_size(config["batch_size"], config["bin_mem"], input[0]),
    log: "logs/kmerBin/{barcode}/shuffle_batch.log"
    benchmark: "benchmarks/kmerBin/{barcode}/shuffle_batch.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        mkdir -p {output}
        sed 1d {input} | shuf > {output}/shuffled
        split -l {params.batch_size} -a3 -d --additional-suffix='.tsv' {output}/shuffled {output}/batch         
        for i in {output}/batch*; do sed -e '1R {input}' -e '1d' -i $i; done
        rm -f {output}/header {output}/shuffled
        """

# kmer binning
#rule umap:
#    input: rules.kmer_freqs.output
#    output: 
#        cluster="kmerBin/{barcode}/hdbscan.tsv",
#	    plot="kmerBin/{barcode}/hdbscan.png",
#    conda: "../envs/kmerBin.yaml"
#    params:
#        n_neighbors = config["umap"]["n_neighbors"],
#        min_dist = config["umap"]["min_dist"],
#        metric = config["umap"]["metric"],
#        n_components = config["umap"]["n_components"],
#	    min_bin_size = config["hdbscan"]["min_bin_size"],
#        min_samples = config["hdbscan"]["min_samples"],
#	    epsilon = config["hdbscan"]["epsilon"],
#    log: "logs/kmerBin/{barcode}/umap.log"
#    benchmark: "benchmarks/kmerBin/{barcode}/umap.txt"
#    threads: config["threads"]["large"]
#    resources:
#        mem_mb = config["bin_mem"] * 1024
#    shell:
#       "NUMBA_NUM_THREADS={threads} python scripts/kmerBin.py -k {input}"
#       " -n {params.n_neighbors} -d {params.min_dist} -r {params.metric} -t {params.n_components}"
#       " -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon}"
#       " -c {output.cluster} -p"
#       " > {log} 2>&1" 

# batch binning
rule umap:
    input: 
        "kmerBin/{barcode}/batch_check/{batch}.tsv",
    output: 
        cluster="kmerBin/{barcode}/{batch}/hdbscan.tsv",
	    plot="kmerBin/{barcode}/{batch}/hdbscan.png",
    conda: "../envs/kmerBin.yaml"
    params:
        n_neighbors = config["umap"]["n_neighbors"],
        min_dist = config["umap"]["min_dist"],
        metric = config["umap"]["metric"],
        n_components = config["umap"]["n_components"],
	    min_bin_size = config["hdbscan"]["min_bin_size"],
        min_samples = config["hdbscan"]["min_samples"],
	    epsilon = config["hdbscan"]["epsilon"],
    log: 
        "logs/kmerBin/{barcode}/{batch}/umap.log"
    benchmark: 
        "benchmarks/kmerBin/{barcode}/{batch}/umap.txt"
    threads: config["threads"]["large"]
    resources:
        mem_mb = config["bin_mem"] * 1024
    shell:
       "NUMBA_NUM_THREADS={threads} python scripts/kmerBin.py -k {input}"
       " -n {params.n_neighbors} -d {params.min_dist} -r {params.metric} -t {params.n_components}"
       " -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

def get_kmerbatch(wildcards):
    batches = glob_wildcards(checkpoints.shuffle_batch.get(**wildcards).output[0] + "/{batch}.tsv").batch
    return expand("kmerBin/{{barcode}}/{batch}/hdbscan.tsv", batch=batches)

rule col_kmerbatch:
    input: get_kmerbatch
    output: "kmerBin/{barcode}/hdbscan.tsv"
    run:
        with open(output[0], "w") as out:
            out.write("read\tlength\tD1\tD2\tbin_id\tbatch_id\n")
            for i in input:
                batch = i.split("/")[-2].replace('batch', 'b')
                with open(i) as f:
                    next(f)
                    for line in f:
                        out.write(line.rstrip() + "\t" + batch + "\n") 

def get_bin(wildcards, pool = True):
    check_val("pool", pool, bool)
    if pool == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_qced(wildcards)
    return expand("kmerBin/{barcode}/hdbscan.tsv", barcode=barcodes)

# split reads by kmerbin
checkpoint cls_kmerbin:
    input: lambda wildcards: get_bin(wildcards, pool = config["pool"])
    output: directory("kmerBin/clusters"),
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            barcode = i.split("/")[-2]
            df_i = pd.read_csv(i, sep="\t")
            # unclustered reads are assigned to bin -1, drop
            df_i = df_i[df_i.bin_id != -1]

            if 'batch_id' in df_i.columns:
                 # concatanate the batch_id bin_id column
                 df_i['bin_id'] = df_i['bin_id'].astype(str) + df_i['batch_id'].astype(str)
            df_i = df_i[["read", "bin_id"]]
            for clust_id, df_clust in df_i.groupby('bin_id'):
                df_clust['read'].to_csv(output[0] + "/{barcode}_{c}.csv".format(barcode=barcode, c=clust_id),
                    header = False, index = False)
        
rule split_bin:
    input: 
        cluster = "kmerBin/clusters/{barcode}_{c}.csv",
        fqs = "qc/qfilt/{barcode}.fastq",
    output: "kmerBin/{barcode}/split/{c}.fastq",
    conda: "../envs/seqkit.yaml"
    log: "logs/kmerBin/{barcode}/split/{c}.log"
    benchmark: "benchmarks/kmerBin/{barcode}/split/{c}.txt"
    shell:
        "seqkit grep {input.fqs} -f {input.cluster} -o {output} 2> {log}"

rule skip_bin:
    input: "qc/qfilt/{barcode}.fastq"
    output: "kmerBin/{barcode}/all.fastq"
    log: "logs/kmerBin/{barcode}/skip_bin.log"
    shell: "cp -p {input} {output} 2> {log}"

# get {barcode} {c} from chekckpoint
def get_kmerBin(wildcards, pool = True, kmerbin = True):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
    
    if kmerbin is True:
        fqs = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            fqs.append("kmerBin/{bc}/split/{kb}.fastq".format(bc=bc, kb=kb))
    else:
        if pool is True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        fqs = expand("kmerBin/{bc}/all.fastq", bc=bcs)
    return fqs
