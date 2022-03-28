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

# barch mode
# shuffle the kmer freqs in batches
checkpoint shuffle_batch:
    input: rules.kmer_freqs.output
    output: directory(OUTPUT_DIR + "/kmerBin/{barcode}/batch")
    conda: "../envs/coreutils.yaml"
    params:
        batch_size = config["batch_size"],
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/shuffle_batch.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/shuffle_batch.txt"
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
rule umap:
    input: rules.kmer_freqs.output
    output: 
        cluster=OUTPUT_DIR + "/kmerBin/{barcode}/hdbscan.tsv",
	    plot=OUTPUT_DIR + "/kmerBin/{barcode}/hdbscan.png",
    conda: "../envs/kmerBin.yaml"
    params:
        n_neighbors = config["umap"]["n_neighbors"],
        min_dist = config["umap"]["min_dist"],
        metric = config["umap"]["metric"],
        n_components = config["umap"]["n_components"],
	    min_bin_size = config["hdbscan"]["min_bin_size"],
        min_samples = config["hdbscan"]["min_samples"],
	    epsilon = config["hdbscan"]["epsilon"],
    log: OUTPUT_DIR + "/logs/kmerBin/{barcode}/umap.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/umap.txt"
    threads: config["threads"]["large"]
    shell:
       "NUMBA_NUM_THREADS={threads} python scripts/kmerBin.py -k {input}"
       " -n {params.n_neighbors} -d {params.min_dist} -r {params.metric} -t {params.n_components}"
       " -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

# batch binning
use rule umap as umap_batch with:
    input: 
        OUTPUT_DIR + "/kmerBin/{barcode}/batch/{batch}.tsv",
    output: 
        cluster=OUTPUT_DIR + "/kmerBin/{barcode}/{batch}/hdbscan.tsv",
	    plot=OUTPUT_DIR + "/kmerBin/{barcode}/{batch}/hdbscan.png",
    log: 
        OUTPUT_DIR + "/logs/kmerBin/{barcode}/{batch}/umap.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/kmerBin/{barcode}/{batch}/umap.txt"
 
def get_kmerbatch(wildcards):
    batches = glob_wildcards(checkpoints.shuffle_batch.get(**wildcards).output[0] + "/{batch}.tsv").batch
    return expand(OUTPUT_DIR + "/kmerBin/{{barcode}}/{batch}/hdbscan.tsv", batch=batches)

rule col_kmerbatch:
    input: get_kmerbatch
    output: OUTPUT_DIR + "/kmerBin/{barcode}/hdbscan_inbatch.tsv"
    run:
        with open(output[0], "w") as out:
            out.write("read\tlength\tD1\tD2\tbin_id\tbatch_id\n")
            for i in input:
                batch = i.split("/")[-2]
                with open(i) as f:
                    next(f)
                    for line in f:
                        out.write(line.rstrip() + "\t" + batch + "\n") 

def get_bin(wildcards, pooling = True, batch_size = -1):
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_qced(wildcards)
    
    out = ""
    if batch_size > 0:
        out = "_inbatch"
    return expand(OUTPUT_DIR + "/kmerBin/{barcode}/hdbscan" + out + ".tsv", barcode=barcodes)

# split reads by kmerbin
checkpoint cls_kmerbin:
    input: lambda wildcards: get_bin(wildcards, pooling = config["pooling"], batch_size = config["batch_size"])
    output: directory(OUTPUT_DIR + "/kmerBin/clusters"),
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
