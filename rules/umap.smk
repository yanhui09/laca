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

# trim UMIs 
rule trim_umi:
    input: rules.nanofilt.output,
    output: OUTPUT_DIR + "/umap/{barcode}/trimmed.fastq",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/trim_umi.log"
    threads: config["threads"]["normal"]
    params:
        f=f5_patterns,
        max_err=config["max_err"],
        min_overlap=config["min_overlap"],
    conda: "../envs/cutadapt.yaml"
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

# kmer calculation
rule kmer_freqs:
    input: rules.trim_umi.output
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
        --ignoreGenomeQuality --SkipMash \
        --clusterAlg single >{log} 2>&1
        """

# racon

# medaka

# call variants

# combine sample-specific amplicons and mapping