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

rule combine_fastq:
    input: INPUT_DIR
    output: temp(OUTPUT_DIR + "/raw_fq/pooled.fastq")
    shell:
        "cat {input}/*.fastq > {output}"

# trim primers 
rule trim_primers:
    input: OUTPUT_DIR + "/raw_fq/{barcode}.fastq"
    output: OUTPUT_DIR + "/umap/{barcode}/trimmed.fastq",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/trim_primers.log"
    threads: config["threads"]["large"]
    params:
        f = f5_patterns,
        max_err = config["max_err"],
        min_overlap = config["min_overlap"],
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
        scripts = "scripts",
        kmer_size = config["kmer_size"],
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
        n_neighbors = config["umap"]["n_neighbors"],
        min_dist = config["umap"]["min_dist"],
        n_components = config["umap"]["n_components"],
	    min_cluster_size = config["hdbscan"]["min_cluster_size"],
        min_samples = config["hdbscan"]["min_samples"],
	    epsilon = config["hdbscan"]["epsilon"],
    threads: config["threads"]["large"]
    shell:
       "NUMBA_NUM_THREADS={threads} python scripts/umap_cluster.py -k {input}"
       " -n {params.n_neighbors} -d {params.min_dist} -t {params.n_components}"
       " -s {params.min_cluster_size} -m {params.min_samples} -e {params.epsilon}"
       " -c {output.cluster} -p"
       " > {log} 2>&1" 

# split reads by cluster
checkpoint cluster_info:
    input: rules.umap.output.cluster,
    output: directory(OUTPUT_DIR + "/umap/{barcode}/clusters"),
    log: OUTPUT_DIR + "/logs/umap/{barcode}/clusters.log"
    params:
        scripts = "scripts"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        df = df[["read", "bin_id"]]
        clusters = df.bin_id.max()
        os.makedirs(output[0])
        for cluster in range(0, clusters+1):
            df.loc[df.bin_id == cluster, "read"].to_csv(output[0] + "/c" + str(cluster) + ".txt", sep="\t", header=False, index=False)
        
rule split_by_cluster:
    input: 
        clusters = OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.txt",
        fastq = rules.nanofilt_umap.output,
    output: OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/clusters/{c}.log"
    conda: "../envs/seqkit.yaml"
    shell:
        "seqkit grep {input.fastq} -f {input.clusters} -o {output} 2> {log}"

# read correction
rule correct_read:
    input:
        fastq = OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    output:
        fasta = OUTPUT_DIR + "/umap/{barcode}/canu_corrected/{c}/{c}.correctedReads.fasta.gz",
    log: OUTPUT_DIR + "/logs/umap/{barcode}/canu/{c}.log"
    threads: config["threads"]["large"]
    conda: "../envs/canu.yaml"
    params:
        amp_size = config["amp_size"],
        prefix = OUTPUT_DIR + "/umap/{barcode}/canu_corrected/{c}",
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
    seqs_folder = checkpoints.split_seqs.get(**wildcards).output[0]
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

# align merged assemblies with raw reads
# reused in racon iterations
rule minimap:
    input: 
      ref = OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/{assembly}.fna",
      fastq = OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    output: OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/{assembly}.paf",
    message: "{wildcards.c} for polish: alignments against {wildcards.assembly} assembly [{wildcards.barcode}]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/polish/minimap/{c}/{assembly}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/{barcode}/polish/minimap/{c}/{assembly}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/raw"
        return(prefix + ".paf", prefix + ".fna")
    else:
        prefix = OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/racon_{iter}".format(barcode=wildcards.barcode, c=wildcards.c, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
        get_racon_input,
    output: OUTPUT_DIR + "/umap/{barcode}/polish/draft/{c}/racon_{iter}.fna"
    message: "Polish {wildcards.c} draft with racon, round={wildcards.iter} [{wildcards.barcode}]"
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR +"/logs/umap/{barcode}/polish/racon/{c}/round{iter}.log"
    benchmark: OUTPUT_DIR +"/benchmarks/umap/{barcode}/polish/racon/{c}/round{iter}.txt"
    threads: 1
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

checkpoint medaka_consensus:
    input:
        fna = expand(OUTPUT_DIR + "/umap/{{barcode}}/polish/draft/{{c}}/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = OUTPUT_DIR + "/umap/{barcode}/clusters/{c}.fastq",
    output: 
        fna = OUTPUT_DIR + "/umap/{barcode}/denoised_seqs/{c}.fna",
    message: "Generate consensus for {wildcards.c} with medaka [{wildcards.barcode}]"
    params:
        m = config["medaka"]["m"],
        _dir = OUTPUT_DIR + "/umap/{barcode}/polish/medaka/{c}",
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/umap/{barcode}/polish/medaka/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/{barcode}/polish/medaka/{c}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fna} -o {params._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        cp {params._dir}/consensus.fasta {output.fna}
        """

# get {barcode} {c} from chekckpoint
def get_denoised_seqs(wildcards, pooling = True):
    if pooling:
        barcodes = ["pooled"]
    elif pooling is False:
        barcodes = glob_wildcards(checkpoints.demultiplex.get(**wildcards).output[0]
        + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    else:
        raise ValueError('Pooling only allows bool type [True/False].\n{} is used in the config file'.format(x))

    fnas = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            fnas.append(OUTPUT_DIR + "/umap/{barcode}/denoised_seqs/{c}.fna".format(barcode=i, c=j))
    return fnas

rule collect_denoised_seqs:
    input: lambda wc: get_denoised_seqs(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/umap/denoised_seqs.fna",
    message: "Collect denoised sequences"
    log: OUTPUT_DIR + "/logs/umap/denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/denoised_seqs.txt"
    shell:
        "cat {input} > {output} 2> {log}"

# to do: call variants with medaka

# dereplicate denoised sequences with mmseqs
rule dereplicate_denoised_seqs:
    input: rules.collect_denoised_seqs.output,
    output: OUTPUT_DIR + "/rep_seqs.fna",
    message: "Dereplicate denoised sequences"
    params:
        tmp = OUTPUT_DIR + "/tmp",
    conda: "../envs/mmseqs2.yaml"
    log: OUTPUT_DIR + "/logs/dereplicate_denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/dereplicate_denoised_seqs.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input} {output} {params.tmp} --threads {threads} --min-seq-id 1 -c 1 > {log} 2>&1"

# create abundance matrix with minimap
rule index:
    input: rules.dereplicate_denoised_seqs.output,
    output: temp(OUTPUT_DIR + "/rep_seqs.mmi")
    message: "Index denoised sequences [Generate abundance matrix]"
    params:
        index_size = config["minimap"]["index_size"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/index.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/index.txt"
    shell: "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule dict:
    input: rules.dereplicate_denoised_seqs.output,
    output: temp(OUTPUT_DIR + "/rep_seqs.dict")
    message: "Dict denoised sequences [Generate abundance matrix]"
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/dict.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/dict.txt"
    shell: "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap_rep_seqs:
    input:
        fq = rules.trim_primers.output,
        mmi = rules.index.output,
        dict = rules.dict.output,
    output: temp(OUTPUT_DIR + "/mapped/{barcode}.bam")
    message: "Re-map trimmed fastq files [Generate abundance matrix]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/minimap/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/minimap/{barcode}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        minimap2 -t {threads} -ax {params.x} {input.mmi} {input.fq} | \
        grep -v "^@" | cat {input.dict} - | \
        samtools view -F 3584 -b - > {output} 2>{log}
        """

rule sort:
    input: rules.minimap_rep_seqs.output
    output: temp(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam")
    params:
        prefix = OUTPUT_DIR + "/mapped/tmp.{barcode}",
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/sort/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/sort/{barcode}.txt"
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m {params.m} -o {output} 2>{log}"

rule samtools_index:        
    input: rules.sort.output
    output: temp(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam.bai")
    params:
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/index/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/index/{barcode}.txt"
    shell:
        "samtools index -m {params.m} -@ 1 {input} {output} 2>{log}"

def get_barcodes(wildcards, type_o):
    barcodes = glob_wildcards(checkpoints.demultiplex.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    if type_o == "bam":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.count", barcode=barcodes)
    else:
        output = barcodes
    return output

# biom format header
rule header_sample:
    input:
        bai = lambda wildcards: get_barcodes(wildcards, "bai"),
    output: temp(OUTPUT_DIR + "/header_sample")
    run:
        with open(output[0], 'w') as f:
            f.write('#OTU ID\t'+ '\t'.join(SAMPLE) + '\n')

rule rowname_seqs:
    input:
        bam = lambda wildcards: get_barcodes(wildcards, "bam"),
        bai = lambda wildcards: get_barcodes(wildcards, "bai"),
    output: temp(OUTPUT_DIR + "/rowname_seqs")
    conda: "../envs/polish.yaml"
    shell:
        """
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 > {output}
        """
       
rule seqs_count:
    input:
        bam = OUTPUT_DIR + "/mapped/{sample}.sorted.bam",
        bai = OUTPUT_DIR + "/mapped/{sample}.sorted.bam.bai"
    output: temp(OUTPUT_DIR + "/mapped/{sample}.count")
    conda: "../envs/polish.yaml"
    shell:
        """
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 > {output}
        """

rule count_matrix:
    input:
        seqs_count = lambda wildcards: get_barcodes(wildcards, "count"),
        header_sample = rules.header_sample.output,
        rowname_seqs = rules.rowname_seqs.output,
    output: OUTPUT_DIR + "/count_matrix.tsv"
    shell:
        """
        paste {input.rowname_seqs} {input.seqs_count} | cat {input.header_sample} - > {output}
        """