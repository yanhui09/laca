# linker and primer info
flinker = config["flinker"]
fprimers = config["fprimer"]
rlinker = config["rlinker"]
rprimers = config["rprimer"]

# reverse complementation
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
flinkerR = revcomp(flinker)
rlinkerR = revcomp(rlinker)
fprimersR = {k: revcomp(v) for k, v in fprimers.items()}
rprimersR = {k: revcomp(v) for k, v in rprimers.items()}

# pattern search for umi using cutadapt
# nanopore possibly sequences either strand
def seqs_join(linker, primer, reverse=False):
    joined = '-g ' + linker + '...' + primer
    if reverse:
        joined = '-G ' + primer + '...' + linker
    return joined
def linked_pattern(linker, primers,reverse=False):
    linked = {k: seqs_join(linker, v, reverse) for k, v in primers.items()}        
    return ' '.join(v for v in linked.values())

# forward
f_pattern1 = linked_pattern(flinker, fprimers)
f_pattern2 = linked_pattern(rlinker, rprimers)
f_pattern = f_pattern1 + ' ' + f_pattern2
# reverse
r_pattern1 = linked_pattern(flinkerR, fprimersR, reverse=True)
r_pattern2 = linked_pattern(rlinkerR, rprimersR, reverse=True)
r_pattern = r_pattern1 + ' ' + r_pattern2
#---------

# trim umi region
rule umi_loc:
    input: OUTPUT_DIR + "/filt_fq/{barcode}.fastq"
    output:
        start=OUTPUT_DIR + "/umi/{barcode}/start.fastq",
        end=OUTPUT_DIR + "/umi/{barcode}/end.fastq",
    log: OUTPUT_DIR + "/logs/umi/umi_loc/{barcode}.log"
    conda: "../envs/seqkit.yaml"
    threads: 1
    params:
        umi_loc=config["umi_loc"]
    shell:
        """
        seqkit subseq -r 1:{params.umi_loc} {input} 2> {log} 1> {output.start}
        seqkit subseq -r -{params.umi_loc}:-1 {input} 2>> {log} 1> {output.end}
        """

# extract UMI sequences 
rule extract_umi:
    input:
        start=rules.umi_loc.output.start,
        end=rules.umi_loc.output.end,
    output:
        umi1=OUTPUT_DIR + "/umi/{barcode}/umi1.fastq",
        umi2=OUTPUT_DIR + "/umi/{barcode}/umi2.fastq",
    log: OUTPUT_DIR + "/logs/umi/extract_umi/{barcode}.log"
    threads: config["threads"]["normal"]
    params:
        f=f_pattern,
        r=r_pattern,
        max_err=config["max_err"],
        min_overlap=config["min_overlap"],
        min_len=config["umi_len"] - config["umi_flexible"],
        max_len=config["umi_len"] + config["umi_flexible"],
    conda: "../envs/cutadapt.yaml"
    shell:
        """
        cutadapt \
            -j {threads} -e {params.max_err} -O {params.min_overlap} \
            -m {params.min_len} -M {params.max_len} \
            --discard-untrimmed \
            {params.f} \
            {params.r} \
            -o {output.umi1} -p {output.umi2} \
            {input.start} {input.end} \
            > {log} 2>&1
        """
        #sleep 10 # some time for writing

# combine UMI sequences
rule concat_umi:
    input:
        umi1=rules.extract_umi.output.umi1,
        umi2=rules.extract_umi.output.umi2,
    output: OUTPUT_DIR + "/umi/{barcode}/umi.fasta"
    log: OUTPUT_DIR + "/logs/umi/concat_umi/{barcode}.log"
    threads: 1
    conda: "../envs/seqkit.yaml"
    shell:
        "seqkit concat {input.umi1} {input.umi2} 2> {log} | seqkit fq2fa -o {output} 2>> {log}"

# cluster UMIs
rule cluster_umi:
    input: rules.concat_umi.output
    output:
        centroid=OUTPUT_DIR + "/umi/{barcode}/centroid.fasta",
        consensus=OUTPUT_DIR + "/umi/{barcode}/consensus.fasta",
        uc=OUTPUT_DIR + "/umi/{barcode}/uc.txt"  
    log: OUTPUT_DIR + "/logs/umi/cluster_umi/{barcode}.log"
    threads: config["threads"]["normal"]
    conda: "../envs/vsearch.yaml"
    params:
        #umi= OUTPUT_DIR + "/umi/{barcode}/umi",
        id=config["umi_id"],
        min_len=2*(config["umi_len"] - config["umi_flexible"]),
        max_len=2*(config["umi_len"] + config["umi_flexible"]),
    shell:
        """
        vsearch \
            --cluster_fast {input} --clusterout_sort -id {params.id} \
            --clusterout_id --sizeout -uc {output.uc}\
            --centroids {output.centroid} --consout {output.consensus} \
            --minseqlength {params.min_len} --maxseqlength {params.max_len} \
            --qmask none --threads {threads} \
            --strand both \
            > {log} 2>&1
        """