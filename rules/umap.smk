# kmer calculation
rule kmer_freqs:
    input: rules.nanofilt.output
    output: OUTPUT_DIR + "/umap/{barcode}/kmer_freqs.txt"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/kmer_freqs.log"
    threads: config["threads"]["normal"]
    conda: "../envs/kmer_freqs.yaml"
    params: 
        scripts="scripts",
	kmer_size=config["kmer_size"],
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
    threads: 1
    conda: "../envs/umap_cluster.yaml"
    params:
        scripts="scripts",
	    size=config["hdbscan_size"],
	    epsilon=config["hdbscan_epsilon"],
    shell:
       "python {params.scripts}/umap_cluster.py -k {input}"
       " -s {params.size} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

# split_by_cluster
rule split_by_cluster:
    input: 
        clusters=rules.umap.output.cluster,
        fastq=rules.nanofilt.output,
    output: OUTPUT_DIR + "/umap/{barcode}/count.txt",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/clusters.log"
    threads: 1
    conda: "../envs/seqkit.yaml"
    params:
        outdir=OUTPUT_DIR + "/umap/{barcode}/clusters",
        scripts="scripts"
    shell:
        """
        bash {params.scripts}/split_by_cluster.sh {input.clusters} {input.fastq} {params.outdir} {output} 2> {log}
        """

# read correction
rule correct_read:
    input:
        counts=rules.split_by_cluster.output,
        fastq=OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    output:
        fasta=OUTPUT_DIR + "/umap/{barcode}/canu_corrected/{c}/{c}.correctedReads.fasta.gz",
        #counts=OUTPUT_DIR + "/umap/{barcode}/counts_canu.txt",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/canu/{c}.log"
    threads: config["threads"]["large"]
    conda: "../envs/canu.yaml"
    params:
        amp_size=config["amp_size"],
        prefix=OUTPUT_DIR + "/umap/{barcode}/canu_corrected/{c}",
    shell:
        """
        canu -correct -p {wildcards.c} -d {params.prefix} \
        -raw -nanopore {input.fastq} \
        genomeSize={params.amp_size} stopOnLowCoverage=1 \
        minInputCoverage=2 minReadLength=500 minOverlapLength=200 \
        useGrid=false maxThreads={threads} > {log} 2>&1
        """

# select draft sequences
checkpoint split_seqs:
    input: rules.correct_read.output.fasta
    output: directory(OUTPUT_DIR + "/umap/{barcode}/canu_corrected/{c}/seqs")
    log: OUTPUT_DIR + "/logs/umap/{barcode}/split_seqs/{c}.log"
    threads: 1
    conda: "../envs/seqkit.yaml"
    shell:
        "seqkit split {input} -i -O {output} > {log} 2>&1"

def get_seqs(wildcards):
    seqs_folder= checkpoints.split_seqs.get(**wildcards).output[0]
    seqs = glob_wildcards(os.path.join(seqs_folder,"{seqs}.fasta.gz")).seqs
    return expand(OUTPUT_DIR + "/umap/{{barcode}}/canu_corrected/{{c}}/seqs/{seqs}.fasta.gz", seqs=seqs)
    
rule drep:
    input: get_seqs
    output: directory(OUTPUT_DIR + "/umap/{barcode}/drep/{c}")
    log: OUTPUT_DIR + "/logs/umap/{barcode}/drep/{c}.log"
    threads: config["threads"]["large"]
    conda: "../envs/drep.yaml"
    params:
        #pa=config["pa"],
        sa=config["sa"],
        S_algorithm=config["S_algorithm"],
        nc=config["nc"],
        minl=config["minl"],
        clusterAlg=config["clusterAlg"],
    shell:
        """
        dRep dereplicate {output} -g {input} \
        -sa {params.sa} \
        --S_algorithm {params.S_algorithm} \
        --clusterAlg {params.clusterAlg} \
        -nc {params.nc} -l {params.minl} -N50W 0 -sizeW 1 \
        --ignoreGenomeQuality --SkipMash >{log} 2>&1
        """

# racon

# medaka

# call variants

# combine sample-specific amplicons and mapping