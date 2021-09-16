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
        fastq=OUTPUT_DIR + "/umap/{barcode}/canu_corrected/{c}.fastq",
        counts=OUTPUT_DIR + "/umap/{barcode}/counts_canu.txt",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/canu/{c}.log"
    threads: config["threads"]["normal"]
    conda: "../envs/canu.yaml"
    params:
        amp_size=config["amp_size"]
        prefix=OUTPUT_DIR + "umap/{barcode}/canu_correctected/{c}"
    shell:
        """
        canu -correct -p {params.prefix}/corrected -raw -nanopore{input.fastq} \
        genomeSize={params.amp_size} stopOnLowCoverage=1 \
        minInputCoverage=2 minReadLength=500 minOverlapLength=200 \

        gunzip {params.prefix}/corrected.fasta.gz
        """
        