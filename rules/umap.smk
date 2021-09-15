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
