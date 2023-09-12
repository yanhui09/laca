check_list_ele("cluster", config["cluster"], ["isONclust", "umap", "isONcorrect", "meshclust"])
localrules: cls_isONclust, cls_umap, cls_meshclust, fqs_split_isONclust, fqs_split_isONclust2, fqs_split_umap, fqs_split_meshclust, prepare_umap_fqs, prepare_meshclust_fqs, fq2fa4meshclust
# "meshclust" is required
if "meshclust" not in config["cluster"]:
    raise ValueError("meshclust is required for cluster")

# isONclust
rule isONclust:
    input: "qc/qfilt/{barcode}.fastq"
    output:
        _dir = temp(directory("clust/isONclust/{barcode}")),
        tsv = temp("clust/isONclust/{barcode}.tsv"),
    params:
        k = config["isONclust"]["k"],
        w = config["isONclust"]["w"],
    conda: "../envs/isONcorCon.yaml"
    log: "logs/clust/isONclust/{barcode}.log"
    benchmark: "benchmarks/clust/isONclust/{barcode}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell:
        """
        # make dummies if failed
        isONclust --k {params.k} --w {params.w} --fastq {input} --outfolder {output._dir} --t {threads} > {log} 2>&1 || true
        if [ -f {output._dir}/final_clusters.tsv ]; then
            mv {output._dir}/final_clusters.tsv {output.tsv}
        else
            mkdir -p {output._dir}
            touch {output.tsv}
        fi
        """

def get_isONclust(wildcards, pool = config["pool"]):
    check_val("pool", pool, bool)
        
    if pool == True:
       bcs = ["pooled"]
    else:
       bcs = get_qced_barcodes(wildcards)
    clusters = expand("clust/isONclust/{bc}.tsv", bc=bcs)
    return clusters

checkpoint cls_isONclust:
    input:
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        clss = lambda wc: get_isONclust(wc),
    output: directory("clust/isONclust/read2cluster")
    params:
        min_size = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.clss):
            # if empty, skip
            if os.stat(i).st_size == 0:
                continue
            barcode = i.split('/')[-1].removesuffix(".tsv")
            df_i = pd.read_csv(i, sep = '\t', header = None)
            df_i.columns = ['cluster', 'seqid']
            for clust_id, df_clust in df_i.groupby('cluster'):
                if len(df_clust) >= params.min_size:
                    df_clust['seqid'].to_csv(output[0] + "/{barcode}_{c1}.csv".format(barcode=barcode, c1=clust_id),
                     header = False, index = False)

rule fqs_split_isONclust:
    input: 
        cluster = "clust/isONclust/read2cluster/{barcode}_{c1}.csv",
        fqs = "qc/qfilt/{barcode}.fastq",
    output: temp("clust/isONclust/split0/{barcode}_{c1}.fastq")
    shell: "seqkit grep -f {input.cluster} {input.fqs} -w0 -o {output} --quiet"

rule fqs_split_isONclust2: 
    input: rules.fqs_split_isONclust.output
    output: temp("clust/isONclust/split/{barcode}_{c1}_0.fastq")
    shell: "cp {input} {output}"

# pseduo fastq if isONclust not used
rule prepare_umap_fqs:
    input: "qc/qfilt/{barcode}.fastq"
    output: temp("clust/umapclust/fqs/{barcode}_0.fastq")
    shell: "cp {input} {output}"

def get_fqs4umap(cluster = config["cluster"]):
    if "isONclust" in cluster:
        fqs_dir = "clust/isONclust/split0"
    else:
        fqs_dir = "clust/umapclust/fqs" 
    return fqs_dir + "/{barcode}_{c1}.fastq"

rule kmer_freqs:
    input: get_fqs4umap()
    output: temp("clust/umapclust/kmer_freqs/{barcode}_{c1}.tsv")
    conda: "../envs/kmerBin.yaml"
    params: 
        kmer_size = config["kmer_size"],
    log: "logs/clust/umapclust/kmer_freqs/{barcode}_{c1}.log"
    benchmark: "benchmarks/clust/umapclust/kmer_freqs/{barcode}_{c1}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    shell:
        "python {workflow.basedir}/scripts/kmer_freqs.py"
        " -k {params.kmer_size}"
        " -r {input} -t {threads}"
        " 2> {log} > {output}"

def get_batch_size(kmer_freq, batch_size = config["batch_size"], mem = config["mem"]["large"]):
    if not os.path.exists(kmer_freq):
        return -1
    # total input read size
    with open(kmer_freq) as f:
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
            est_size = int((1024 * mem - 254) / m)
            if est_size > total:
                batch_size = total
            else:
                batch_size = est_size
    return batch_size

rule umapclust:
    input: rules.kmer_freqs.output
    output: 
        batch = temp(directory("clust/umapclust/{barcode}_{c1}")),
        read2cluster = "clust/umapclust/{barcode}_{c1}.tsv",
    conda: "../envs/kmerBin.yaml"
    params:
        bc_clss= "{barcode}_{c1}",
        batch = lambda wc, input: get_batch_size(kmer_freq = input[0]),
        n_neighbors = config["umap"]["n_neighbors"],
        min_dist = config["umap"]["min_dist"],
        metric = config["umap"]["metric"],
        n_components = config["umap"]["n_components"],
	    min_bin_size = config["hdbscan"]["min_bin_size"],
        min_samples = config["hdbscan"]["min_samples"],
	    epsilon = config["hdbscan"]["epsilon"],
    log: "logs/clust/umapclust/{barcode}_{c1}/umapclust.log",
    benchmark: "benchmarks/clust/umapclust/{barcode}_{c1}/umapclust.txt",
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["default"],
    shell:
        """
        mkdir -p {output.batch}
        sed 1d {input} | shuf --random-source=<(yes 123) > {output.batch}/shuffled
        split -l {params.batch} -a3 -d --additional-suffix='.tsv' {output.batch}/shuffled {output.batch}/batch         
        for i in {output.batch}/batch*; do sed -e '1R {input}' -e '1d' $i > /tmp/kb_batch && cat /tmp/kb_batch > $i; done
        rm -f {output.batch}/header {output.batch}/shuffled

        # umap cluster in loop
        for i in {output.batch}/batch*tsv; do
            batchid=$(basename $i .tsv)
            
            NUMBA_NUM_THREADS={threads} python {workflow.basedir}/scripts/umapclust.py -k $i \
            -n {params.n_neighbors} -d {params.min_dist} -r {params.metric} -t {params.n_components} \
            -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon} \
            -c  {output.batch}/umap_$batchid.tsv \
            > {log} 2>&1
            rm -f $i
        done
        
        # combine all batch and add batchid column to read2cluster
        python {workflow.basedir}/scripts/combine_umapclust.py {output.read2cluster} {output.batch}
        """

def get_umapclust(wildcards, pool = config["pool"], cluster = config["cluster"], fqs=False):
    check_val("pool", pool, bool)
    
    clusters = []
    fqs = []

    if "isONclust" in cluster:
        bc_clss = glob_wildcards(checkpoints.cls_isONclust.get(**wildcards).output[0] + "/{bc_cls}.csv").bc_cls
        for i in bc_clss:
            bc, clss = i.split("_")
            clusters.append("clust/umapclust/{bc}_{clss}.tsv".format(bc=bc, clss=clss))
            fqs.append("clust/isONclust/split/{bc}_{clss}.fastq".format(bc=bc, clss=clss))
    else:
        if pool == True:
           bcs = ["pooled"]
        else:
           bcs = get_qced_barcodes(wildcards)
        clusters = expand("clust/umapclust/{bc}_0.tsv", bc=bcs)
        fqs = expand("clust/umapclust/fqs/{bc}_0.fastq", bc=bcs)
    if fqs is True:
        return fqs
    else:
        return clusters

# split reads by umap cluster
checkpoint cls_umap:
    input:
        ".qc_DONE",
        lambda wc: sorted(set(get_filt(wc) + expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)))),
        lambda wc: get_umapclust(wc, fqs=True), 
        clusters = lambda wc: get_umapclust(wc),
    output: directory("clust/umapclust/read2cluster")
    params:
        min_size = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.clusters):
            # if empty, skip
            if os.stat(i).st_size == 0:
                continue
            bc_cls1 = i.split('/')[-1].removesuffix(".tsv")
            df_i = pd.read_csv(i, sep="\t")
            # unclustered reads are assigned to bin -1, drop
            df_i = df_i[df_i.bin_id != -1]

            if 'batch_id' in df_i.columns:
                 # concatanate the batch_id bin_id column
                 df_i['bin_id'] = df_i['bin_id'].astype(str) + df_i['batch_id'].astype(str)
            df_i = df_i[["read", "bin_id"]]
            for cluster_id, df_clust in df_i.groupby('bin_id'):
                if len(df_clust) >= params.min_size:
                    df_clust['read'].to_csv(output[0] + "/{bc_cls1}_{cls2}.csv".format(bc_cls1=bc_cls1, cls2=cluster_id), header = False, index = False)

use rule fqs_split_isONclust as fqs_split_umap with:
    input:
        cluster = "clust/umapclust/read2cluster/{barcode}_{c1}_{c2}.csv",
        fqs = "qc/qfilt/{barcode}.fastq",
    output: 
        temp("clust/umapclust/split/{barcode}_{c1}_{c2}.fastq")

# meshclust
use rule prepare_umap_fqs as prepare_meshclust_fqs with:
    input: 
        "qc/qfilt/{barcode}.fastq"
    output: 
        temp("clust/meshclust/fqs/{barcode}_0_0.fastq")

def get_fqs4isONcorrect(cluster = config["cluster"]):
    if "umap" in cluster:
        fqs_dir = "clust/umapclust/split"
    elif "isONclust" in cluster:
        fqs_dir = "clust/isONclust/split"
    else:
        fqs_dir = "clust/meshclust/fqs"
    return fqs_dir + "/{barcode}_{c1}_{c2}.fastq"

# run_isoncorrect provide multithread processing for isONcorrect in batches
rule isONcorrect:
    input: get_fqs4isONcorrect()
    output:
        rundir = temp(directory("clust/isONcorrect/{barcode}_{c1}/{c2}")), 
        fastq = temp("clust/isONcorrect/{barcode}_{c1}_{c2}.fastq")
    conda: "../envs/isONcorCon.yaml"
    params:
        max_seqs = 2000,
    log: "logs/clust/isONcorrect/{barcode}_{c1}_{c2}.log"
    benchmark: "benchmarks/clust/isONcorrect/{barcode}_{c1}_{c2}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    shell:
        """
        mkdir -p {output.rundir}
        cp {input} {output.rundir}/0.fastq
        # number of reads in fastqs
        nlines=$(cat {output.rundir}/0.fastq | wc -l)
        if [ $nlines -gt $((4*{params.max_seqs})) ]; then
            run_isoncorrect --t {threads} --fastq_folder {output.rundir} --outfolder {output.rundir} --set_w_dynamically --split_wrt_batches --max_seqs {params.max_seqs} > {log} 2>&1
        else
            run_isoncorrect --t {threads} --fastq_folder {output.rundir} --outfolder {output.rundir} --set_w_dynamically --max_seqs {params.max_seqs} > {log} 2>&1
        fi
        mv {output.rundir}/0/corrected_reads.fastq {output.fastq}
        """

def get_fqs4meshclust(cluster = config["cluster"]):
    if "isONcorrect" in cluster:
        return rules.isONcorrect.output.fastq
    else:
        return get_fqs4isONcorrect(cluster)

rule fq2fa4meshclust:
    input: get_fqs4meshclust()
    output: temp("clust/meshclust/{barcode}_{c1}_{c2}.fasta")
    shell: "seqkit fq2fa {input} -w0 -o {output} --quiet"

rule meshclust:
    input: "clust/meshclust/{barcode}_{c1}_{c2}.fasta"
    output: "clust/meshclust/{barcode}_{c1}_{c2}.tsv"
    params:
        t = "-t " + str(config["meshclust"]["t"]) if config["meshclust"]["t"] else "",
    singularity: "docker://yanhui09/identity:latest"
    log: "logs/clust/meshclust/{barcode}_{c1}_{c2}.log"
    benchmark: "benchmarks/mechclust/{barcode}_{c1}_{c2}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell: 
        "meshclust -d {input} -o {output} -c {threads} {params.t} > {log} 2>&1"

def get_meshclust(wildcards, cluster = config["cluster"], fastq = False):
    if "umap" in cluster: 
        bc_clss = glob_wildcards(checkpoints.cls_umap.get(**wildcards).output[0] + "/{bc_clss}.csv").bc_clss
        fqs_dir = "clust/umapclust/split"
    elif "isONclust" in cluster:
        bc_cls1 = glob_wildcards(checkpoints.cls_isONclust.get(**wildcards).output[0] + "/{bc_cls1}.csv").bc_cls1
        bc_clss = [i + "_0" for i in bc_cls1]
        fqs_dir = "clust/isONclust/split"
    else:
        bcs = get_qced_barcodes(wildcards)
        bc_clss = [i + "_0_0" for i in bcs]
        fqs_dir = "clust/meshclust/fqs"
    
    if "isONcorrect" in cluster:
        fqs_dir = "clust/isONcorrect"
    
    if fastq == True:
        return expand(fqs_dir + "/{bc_cls}.fastq", bc_cls=bc_clss)
    else: 
        return expand("clust/meshclust/{bc_cls}.tsv", bc_cls=bc_clss)

checkpoint cls_meshclust:
    input:  
        ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: get_meshclust(wc, fastq = True),
        cluster = lambda wc: get_meshclust(wc, fastq = False),
    output: directory("clust/clusters")
    params:
        min_size = 4,
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.cluster):
            # if empty, skip
            if os.stat(i).st_size == 0:
                continue
            barcode, c1, c2 = [ i.split('/')[-1].removesuffix(".tsv").split('_')[index] for index in [-3, -2, -1] ]
            df_i = pd.read_csv(i, sep = '\t', header = None)
            df_i.columns = ['cluster', 'seqid', 'identity', 'status']
            # trim '>' in seqid
            df_i['seqid'] = df_i['seqid'].str[1:] 
            for clust_id, df_clust in df_i.groupby('cluster'):
                if len(df_clust) >= params.min_size:
                    df_clust['seqid'].to_csv(output[0] + "/{barcode}_{c1}_{c2}_{c3}.csv".format(barcode=barcode, c1=c1, c2=c2, c3=clust_id),
                     header = False, index = False)
                    # extract ref "status = C" to csv file
                    df_clust[df_clust['status'] == 'C']['seqid'].to_csv(output[0] + "/{barcode}_{c1}_{c2}_{c3}.centroid".format(barcode=barcode, c1=c1, c2=c2, c3=clust_id),
                     header = False, index = False)

rule fqs_split_meshclust:
    input:
        members = "clust/clusters/{barcode}_{c1}_{c2}_{c3}.csv",
        centroid = "clust/clusters/{barcode}_{c1}_{c2}_{c3}.centroid",
        split = get_fqs4meshclust(),
    output:
        members = temp("clust/members/{barcode}_{c1}_{c2}_{c3}.fastq"),
        centroid = temp("clust/centroids/{barcode}_{c1}_{c2}_{c3}.fasta"),
    shell:
        """
        seqkit grep -f {input.members} {input.split} -o {output.members} --quiet
        seqkit grep -f {input.centroid} {input.split} --quiet | seqkit fq2fa | \
        seqkit replace -p '^(.+)$' -r '{wildcards.barcode}_{wildcards.c1}_{wildcards.c2}_{wildcards.c3}' -o {output.centroid} --quiet
        """

def get_kmerclust(wildcards, export_centroids = False):

    bc_clss = glob_wildcards(checkpoints.cls_meshclust.get(**wildcards).output[0] + "/{bc_clss}.csv").bc_clss
    members = expand("clust/members/{bc_clss}.fastq", bc_clss=bc_clss)
    centroids = expand("clust/centroids/{bc_clss}.fasta", bc_clss=bc_clss)
    if export_centroids == True:
        return members + centroids
    else:
        return members
