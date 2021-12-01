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

# trim primers 
rule trim_primers:
    input: rules.collect_fastq.output,
    output: OUTPUT_DIR + "/umap/{barcode}/trimmed.fastq",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/trim_primers.log"
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

# quality filter
rule nanofilt_umap:
    input:  rules.trim_primers.output,
    output: OUTPUT_DIR + "/umap/{barcode}/filt.fastq"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/nanofilt.log"
    conda: "../envs/nanofilt.yaml"
    shell: 
        """
        cat {input} | NanoFilt -q 8 -l 800 --maxlength 2000 2> {log} 1> {output}
        """

# kmer calculation
rule kmer_freqs:
    input: rules.nanofilt_umap.output
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
checkpoint split_by_cluster:
    input: 
        clusters=rules.umap.output.cluster,
        fastq=rules.nanofilt_umap.output,
    output: directory(OUTPUT_DIR + "/umap/{barcode}/clusters"),
    log: OUTPUT_DIR + "/logs/umap/{barcode}/clusters.log"
    conda: "../envs/seqkit.yaml"
    params:
        scripts="scripts"
    shell:
        "bash {params.scripts}/split_by_cluster.sh {input.clusters} {input.fastq} {output} 2> {log}"

# read correction
rule correct_read:
    input:
        fastq=OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    output:
        fasta=OUTPUT_DIR + "/umap/{barcode}/canu_corrected/{c}/{c}.correctedReads.fasta.gz",
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
    conda: "../envs/seqkit.yaml"
    shell:
        "seqkit split {input} -i -O {output} > {log} 2>&1"

def get_seqs(wildcards):
    seqs_folder= checkpoints.split_seqs.get(**wildcards).output[0]
    seqs = glob_wildcards(os.path.join(seqs_folder,"{seqs}.fasta.gz")).seqs
    return expand(OUTPUT_DIR + "/umap/{{barcode}}/canu_corrected/{{c}}/seqs/{seqs}.fasta.gz", seqs=seqs)
    
rule fastani:
    input: get_seqs
    output: OUTPUT_DIR + "/umap/{barcode}/fastani/{c}.txt"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/fastani/{c}.log"
    threads: config["threads"]["large"]
    conda: "../envs/fastani.yaml"
    params:
        _dir = OUTPUT_DIR + "/umap/{barcode}/fastani",
        k = config["fastANI"]["k"],
        fragLen = config["fastANI"]["fragLen"],
        minFrac = config["fastANI"]["minFrac"],
    shell:
        """
        ls {input} > {params._dir}/{wildcards.c}.list
        fastANI \
        --ql {params._dir}/{wildcards.c}.list \
        --rl {params._dir}/{wildcards.c}.list \
        --output {output} \
        --threads {threads} \
        --kmer {params.k} \
        --fragLen {params.fragLen} \
        --minFraction {params.minFrac} \
        > {log} 2>&1
        """

rule inter_ANI:
    input: rules.fastani.output,
    output: OUTPUT_DIR + "/umap/{barcode}/fastani/{c}_interANI.txt",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/fastani/{c}_interANI.log"
    run:
        import pandas as pd
        f = pd.read_csv(input[0], sep="\t", header=None)
        f.columns = ["query", "ref", "ANI", "bi_frag", "total_frag"]
        f.drop_duplicates(subset=["query", "ref"], keep=False)
        f.groupby(["query"]).ANI.mean().sort_values(ascending=False).to_csv(output[0], sep="\t", header=False, index=True)

# use the sequence with highest inter-ANI as draft to polish
# rename header to avoid incompatibility with racon
rule choose_draft:
    input: rules.inter_ANI.output,
    output: OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/raw.fna",
    run:
        import pandas as pd
        import gzip
        f = pd.read_csv(input[0], sep="\t", header=None)
        f1 = f.iloc[0,0]
        with open(output[0], "w") as out:
            with gzip.open(f1, "rt") as inp:
                for line in inp:
                    if line.startswith(">"):
                        line = ">" + wildcards.barcode + "_" + wildcards.c + "\n"
                    out.write(line)
