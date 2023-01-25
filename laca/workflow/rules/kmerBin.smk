rule kmer_freqs:
    input: "qc/qfilt/{barcode}.fastq"
    output: temp("kmerBin/{barcode}/kmer_freqs.txt")
    conda: "../envs/kmerBin.yaml"
    params: 
        kmer_size = config["kmer_size"],
    log: "logs/kmerBin/kmer_freqs/{barcode}.log"
    benchmark: "benchmarks/kmerBin/kmer_freqs/{barcode}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    shell:
        "python {workflow.basedir}/scripts/kmerFreqs.py"
        " -k {params.kmer_size}"
        " -r {input} -t {threads}"
        " 2> {log} > {output}"

def get_batch_size(kmer_file, batch_size = config["batch_size"], ram = config["bin_mem"]):
    if not os.path.exists(kmer_file):
        return -1
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
            # estimated constan m
            # max_rss y = 254 + 0.0514x R2=0.96 
            # max_vms y = 1570 + 0.0969x R2 = 0.81
            
            m = 0.06
            est_size = int((1024 * ram - 254) / m)
            if est_size > total:
                batch_size = total
            else:
                batch_size = est_size
    return batch_size

checkpoint shuffle_batch:
    input: rules.kmer_freqs.output
    output: temp(directory("kmerBin/{barcode}/batch_check"))
    conda: "../envs/coreutils.yaml"
    params:
        batch_size = lambda wc, input: get_batch_size(kmer_file = input[0]),
    benchmark: "benchmarks/kmerBin/shuffle_batch/{barcode}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        mkdir -p {output}
        sed 1d {input} | shuf --random-source=<(yes 123) > {output}/shuffled
        split -l {params.batch_size} -a3 -d --additional-suffix='.tsv' {output}/shuffled {output}/batch         
        for i in {output}/batch*; do sed -e '1R {input}' -e '1d' -i $i; done
        rm -f {output}/header {output}/shuffled
        """

rule umap:
    input: 
        "kmerBin/{barcode}/batch_check/{batch}.tsv",
    output: 
        cluster=temp("kmerBin/{barcode}/{batch}/hdbscan.tsv"),
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
        "logs/kmerBin/umap/{barcode}_{batch}.log"
    benchmark: 
        "benchmarks/kmerBin/umap/{barcode}_{batch}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["bin_mem"],
        time = config["runtime"]["default"],
    shell:
       "NUMBA_NUM_THREADS={threads} python {workflow.basedir}/scripts/kmerBin.py -k {input}"
       " -n {params.n_neighbors} -d {params.min_dist} -r {params.metric} -t {params.n_components}"
       " -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

def get_kmerbatch(wildcards):
    batches = glob_wildcards(checkpoints.shuffle_batch.get(**wildcards).output[0] + "/{batch}.tsv").batch
    return expand("kmerBin/{{barcode}}/{batch}/hdbscan.tsv", batch=batches)

localrules: col_kmerbatch, cls_kmerbin, fqs_split
rule col_kmerbatch:
    input: 
        "kmerBin/{barcode}/batch_check",
        batch = lambda wc: get_kmerbatch(wc)
    output: "kmerBin/{barcode}/hdbscan.tsv"
    run:
        with open(output[0], "w") as out:
            out.write("read\tlength\tD1\tD2\tbin_id\tbatch_id\n")
            for i in input.batch:
                batch = i.split("/")[-2].replace('batch', 'b')
                with open(i) as f:
                    next(f)
                    for line in f:
                        out.write(line.rstrip() + "\t" + batch + "\n") 

def get_bin(wildcards, pool = config["pool"]):
    check_val("pool", pool, bool)
    if pool == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_qced_barcodes(wildcards)
    return expand("kmerBin/{barcode}/hdbscan.tsv", barcode=barcodes)

# split reads by kmerbin
checkpoint cls_kmerbin:
    input:
        ".qc_DONE",
        lambda wc: sorted(set(get_filt(wc) + expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)))), 
        bin = lambda wc: get_bin(wc)
    output: directory("kmerBin/clusters"),
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.bin):
            barcode = i.split("/")[-2]
            df_i = pd.read_csv(i, sep="\t")
            # unclustered reads are assigned to bin -1, drop
            df_i = df_i[df_i.bin_id != -1]

            if 'batch_id' in df_i.columns:
                 # concatanate the batch_id bin_id column
                 df_i['bin_id'] = df_i['bin_id'].astype(str) + df_i['batch_id'].astype(str)
            df_i = df_i[["read", "bin_id"]]
            for clust_id, df_clust in df_i.groupby('bin_id'):
                df_clust['read'].to_csv(output[0] + "/{barcode}_{c}.csv".format(
                    barcode=barcode, c=clust_id), header = False, index = False)
rule fqs_split:
    input:
        cluster = "kmerBin/clusters/{barcode}_{c}.csv",
        fqs = "qc/qfilt/{barcode}.fastq",
    output: temp("kmerBin/split/{barcode}_{c}_0.fastq"),
    log: "logs/kmerBin/fqs_split/{barcode}_{c}.log"
    benchmark: "benchmarks/kmerBin/fqs_split/{barcode}_{c}.txt"
    shell: "seqkit grep -f {input.cluster} {input.fqs} -w0 -o {output} --quiet 2> {log}"

# get {barcode} {c} from chekckpoint
def get_kmerBin(wildcards, pool = config["pool"], kmerbin = config["kmerbin"]):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
    
    if kmerbin is True:
        fqs = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            fqs.append("kmerBin/split/{bc}_{kb}_0.fastq".format(bc=bc, kb=kb))
    else:
        if pool is True:
           bcs = ["pooled"]
        else:
           bcs = get_qced_barcodes(wildcards)
        fqs = expand("qc/qfilt/{bc}.fastq", bc=bcs)
    return fqs
