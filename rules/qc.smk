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
    input: rules.collect_fastq.output
    output: OUTPUT_DIR + "/raw/primers_trimmed/{barcode}.fastq"
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5_patterns,
        max_err = config["max_err"],
        min_overlap = config["min_overlap"],
    log: OUTPUT_DIR + "/logs/raw/{barcode}/trim_primers.log"
    benchmark: OUTPUT_DIR + "/benchmarks/raw/{barcode}/trim_primers.txt"
    threads: config["threads"]["large"]
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
rule q_filter:
    input:  rules.trim_primers.output,
    output: OUTPUT_DIR + "/raw/qfilt/{barcode}.fastq"
    conda: "../envs/seqkit.yaml"
    params:
        Q = config["seqkit"]["min-qual"],
        m = config["seqkit"]["min-len"],
        M = config["seqkit"]["max-len"],
    log: OUTPUT_DIR + "/logs/raw/{barcode}/q_filter.log"
    benchmark: OUTPUT_DIR + "/benchmarks/raw/{barcode}/q_filter.txt"
    threads: config["threads"]["normal"]
    shell: "seqkit seq -j {threads} -Q {params.Q} -m {params.m} -M {params.M} {input} > {output} 2> {log}"

#  pooling fqs for sensitivity 
rule combine_fastq:
    input: lambda wc: expand(OUTPUT_DIR + "/raw/qfilt/{barcode}.fastq", barcode=get_demultiplexed(wc))
    output: temp(OUTPUT_DIR + "/raw/qfilt/pooled.fastq")
    shell:
        "cat {input} > {output}"

def get_filt(wildcards, pooling = True):
    barcodes = get_demultiplexed(wildcards) 
    check_val_pool(pooling)
    if pooling == True:
        barcodes.append("pooled")
    return expand(OUTPUT_DIR + "/raw/qfilt/{barcode}.fastq", barcode=barcodes)
