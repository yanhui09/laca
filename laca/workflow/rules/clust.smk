check_list_ele("cluster", config["cluster"], ["isONclust", "umapclust", "meshclust"])
localrules: cls_isONclust, cls_umapclust, cls_meshclust, fqs_split_isONclust, fqs_split_isONclust2, fqs_split_umapclust, fqs_split_meshclust, prepare_umapclust_fqs, prepare_meshclust_fqs, fq2fa4meshclust
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
    output: temp("clust/isONclust/{barcode}_{c1}.split.fastq")
    shell: "seqkit grep -f {input.cluster} {input.fqs} -w0 -o {output} --quiet"

rule fqs_split_isONclust2: 
    input: rules.fqs_split_isONclust.output
    output: temp("clust/isONclust/split/{barcode}_{c1}_0.fastq")
    shell: "cp {input} {output}"

# pseduo fastq if isONclust not used
rule prepare_umapclust_fqs:
    input: "qc/qfilt/{barcode}.fastq"
    output: temp("clust/umapclust/fqs/{barcode}_0.fastq")
    shell: "cp {input} {output}"

def get_fqs4umapclust(cluster = config["cluster"]):
    if "isONclust" in cluster:
        fastq = "clust/isONclust/{barcode}_{c1}.split.fastq"
    else:
        fastq = "clust/umapclust/fqs/{barcode}_{c1}.fastq" 
    return fastq

rule kmer_freqs:
    input: get_fqs4umapclust()
    output: temp("clust/umapclust/{barcode}_{c1}.kmerfreq")
    conda: "../envs/umapclust.yaml"
    params: 
        kmer_size = 5,
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

def get_batch_size(batch_size = config["umapclust"]["max_batch_size"], mem = config["mem"]["large"]):
    try: 
        batch_size = int(batch_size)
        if batch_size <= 0:
            batch_size = -1
    # auto estimation of batch size based on RAM
    except:
        if batch_size == "auto":
            # estimated constan m
            # max_rss y = 254 + 0.0514x R2=0.96 
            # max_vms y = 1570 + 0.0969x R2 = 0.81
            m = 0.06
            est_size = int((1024 * mem - 254) / m)
    return batch_size

rule umapclust:
    input: rules.kmer_freqs.output
    output: "clust/umapclust/{barcode}_{c1}.tsv",
    conda: "../envs/umapclust.yaml"
    params:
        prefix = "clust/umapclust/{barcode}_{c1}",
        max_batch_size = get_batch_size(),
        n_neighbors = config["umapclust"]["umap"]["n_neighbors"],
        min_dist = config["umapclust"]["umap"]["min_dist"],
        metric = "cosine",
        n_components = 2,
	    min_bin_size = config["umapclust"]["hdbscan"]["min_bin_size"],
        min_samples = config["umapclust"]["hdbscan"]["min_samples"],
	    epsilon = config["umapclust"]["hdbscan"]["epsilon"],
    log: "logs/clust/umapclust/{barcode}_{c1}/umapclust.log",
    benchmark: "benchmarks/clust/umapclust/{barcode}_{c1}/umapclust.txt",
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["default"],
    shell:
        """
        # number of reads in fastqs
        nlines=$(cat {input} | wc -l)
        # header excluded
        nlines=$((nlines - 1))
        if [ {params.max_batch_size} -eq -1 ] || [ $nlines -le {params.max_batch_size} ]; then
            # if batch_size <= 0 or batch_size > nlines, run umapclust directly
            NUMBA_NUM_THREADS={threads} python {workflow.basedir}/scripts/umapclust.py -k {input} \
            -n {params.n_neighbors} -d {params.min_dist} -r {params.metric} -t {params.n_components} \
            -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon} \
            -c {output} \
            > {log} 2>&1
        else
            # determine the minimum partion size (number of batches), ceiling division
            min_part_size=$(((nlines + {params.max_batch_size} - 1) / {params.max_batch_size}))
            # determine the number of lines per batch, ceiling division
            nlines_per_batch=$(((nlines + min_part_size - 1) / min_part_size))
            # if batch folder not exist, mkdir, shuffle & split
            if [ ! -d {params.prefix} ]; then
              mkdir -p {params.prefix}
              # shuffle split fastq by barcode if in pooling mode
              sed 1d {input} | shuf --random-source=<(yes 123) > {params.prefix}/shuffled
              split -l $nlines_per_batch -a3 -d --additional-suffix='.tsv' {params.prefix}/shuffled {params.prefix}/b >{log} 2>&1         
              for i in {params.prefix}/b*.tsv; do sed -e '1R {input}' -e '1d' $i > /tmp/kb_batch && cat /tmp/kb_batch > $i; done
              rm -f {params.prefix}/shuffled
            fi
            # umap cluster in loop
            for i in {params.prefix}/b*.tsv; do
                batchid=$(basename $i .tsv)
                
                NUMBA_NUM_THREADS={threads} python {workflow.basedir}/scripts/umapclust.py -k $i \
                -n {params.n_neighbors} -d {params.min_dist} -r {params.metric} -t {params.n_components} \
                -s {params.min_bin_size} -m {params.min_samples} -e {params.epsilon} \
                -c  {params.prefix}/umapclust_$batchid.tsv \
                >> {log} 2>&1
                rm -f $i
            done
            # combine all batch and add batchid column to read2cluster
            python {workflow.basedir}/scripts/combine_umapclust.py {output} {params.prefix}
            rm -rf {params.prefix}
        fi
        """

def get_umapclust(wildcards, pool = config["pool"], cluster = config["cluster"], fqs=False):
    check_val("pool", pool, bool)
    
    clusters = []
    fqs = []

    if "isONclust" in cluster:
        bc_clss = glob_wildcards(checkpoints.cls_isONclust.get(**wildcards).output[0] + "/{bc_cls}.csv").bc_cls
        for i in bc_clss:
            clusters.append("clust/umapclust/{bc_cls}.tsv".format(bc_cls=i))
            # if isONclust + umapclust
            fqs.append("clust/isONclust/{bc_cls}.split.fastq".format(bc_cls=i))
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
checkpoint cls_umapclust:
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

use rule fqs_split_isONclust as fqs_split_umapclust with:
    input:
        cluster = "clust/umapclust/read2cluster/{barcode}_{c1}_{c2}.csv",
        fqs = "qc/qfilt/{barcode}.fastq",
    output: 
        temp("clust/umapclust/split/{barcode}_{c1}_{c2}.fastq")

# meshclust
use rule prepare_umapclust_fqs as prepare_meshclust_fqs with:
    input: 
        "qc/qfilt/{barcode}.fastq"
    output: 
        temp("clust/meshclust/fqs/{barcode}_0_0.fastq")

def get_fqs4meshclust(cluster = config["cluster"]):
    if "umapclust" in cluster:
        fqs_dir = "clust/umapclust/split"
    elif "isONclust" in cluster:
        fqs_dir = "clust/isONclust/split"
    else:
        fqs_dir = "clust/meshclust/fqs"
    return fqs_dir + "/{barcode}_{c1}_{c2}.fastq"

rule fq2fa4meshclust:
    input: get_fqs4meshclust()
    output: temp("clust/meshclust/{barcode}_{c1}_{c2}.fasta")
    shell: "seqkit fq2fa {input} -w0 -o {output} --quiet"

rule meshclust:
    input: "clust/meshclust/{barcode}_{c1}_{c2}.fasta"
    output: "clust/meshclust/{barcode}_{c1}_{c2}.tsv"
    params:
        prefix = "clust/meshclust/{barcode}_{c1}_{c2}",
        t = "-t " + str(config["meshclust"]["t"]) if config["meshclust"]["t"] else "",
        max_batch_size = -1 if int(config["meshclust"]["max_batch_size"]) <= 0 else int(config["meshclust"]["max_batch_size"]) * 2,
    singularity: "docker://yanhui09/identity:latest"
    log: "logs/clust/meshclust/{barcode}_{c1}_{c2}.log"
    benchmark: "benchmarks/mechclust/{barcode}_{c1}_{c2}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell: 
        """
        if [ {params.max_batch_size} -eq -1 ]; then
          meshclust -d {input} -o {output} -c {threads} {params.t} > {log} 2>&1
        else
          # number of reads in fastqs
          nlines=$(cat {input} | wc -l)
          if [ $nlines -le {params.max_batch_size} ]; then
            meshclust -d {input} -o {output} -c {threads} {params.t} > {log} 2>&1
          else
            # determine the minimum partion size (number of batches), ceiling division
            min_part_size=$(((nlines + {params.max_batch_size} - 1) / {params.max_batch_size}))
            # determine the number of lines per batch, ceiling division, the nearest multiples of 2
            nlines_per_batch=$(((nlines + min_part_size * 2 - 1) / (min_part_size * 2) * 2))
            # if batches folder not exist, mkdir & split
            if [ ! -d {params.prefix} ]; then
              mkdir -p {params.prefix}
              split -l $nlines_per_batch -a3 -d --additional-suffix='.fasta' {input} {params.prefix}/b >{log} 2>&1
            fi
            for fa in {params.prefix}/b*.fasta; do
              batch_id=$(basename $fa | cut -d'.' -f1)
              if [ -f {params.prefix}/$batch_id.tsv ]; then
                continue
              fi
              meshclust -d $fa -o {params.prefix}/$batch_id.tsv -c {threads} {params.t} >> {log} 2>&1
              # add batchid after cluster (the first column)
              sed -i "s/\t/$batch_id\t/" {params.prefix}/$batch_id.tsv
              rm -f $fa
            done
            cat {params.prefix}/b*.tsv > {output}
            rm -rf {params.prefix} 
          fi 
        fi
        """

def get_meshclust(wildcards, cluster = config["cluster"], pool = config["pool"]):
    if "umapclust" in cluster: 
        bc_clss = glob_wildcards(checkpoints.cls_umapclust.get(**wildcards).output[0] + "/{bc_clss}.csv").bc_clss
    elif "isONclust" in cluster:
        bc_cls1 = glob_wildcards(checkpoints.cls_isONclust.get(**wildcards).output[0] + "/{bc_cls1}.csv").bc_cls1
        bc_clss = [i + "_0" for i in bc_cls1]
    else:
        if pool == True:
            bcs = ["pooled"]
        else:
            bcs = get_qced_barcodes(wildcards)
        bc_clss = [i + "_0_0" for i in bcs]
    return expand("clust/meshclust/{bc_cls}.tsv", bc_cls=bc_clss)

checkpoint cls_meshclust:
    input:  
        ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        cluster = lambda wc: get_meshclust(wc),
    output: directory("clust/clusters")
    params:
        min_size = config["min_cluster_size"],
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
        split = "qc/qfilt/{barcode}.fastq",
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
