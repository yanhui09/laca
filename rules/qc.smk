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

rule subsample:
    input: rules.collect_fastq.output
    output:
        p = temp(OUTPUT_DIR + "/qc/subsampled/{barcode}_p.fastq"),
        n = temp(OUTPUT_DIR + "/qc/subsampled/{barcode}.fastq"),
    conda: "../envs/seqkit.yaml"
    params:
        p = config["seqkit"]["p"],
        n = config["seqkit"]["n"],
    log: OUTPUT_DIR + "/logs/subsample/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/subsample/{barcode}.txt"
    threads: 1
    shell:
        # I don't know why pipe fails, 
        # portion subsampling (required by seqkit sample by number) as temp file instead.
        """
        seqkit sample -p {params.p} -j {threads} {input} -o {output.p} 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -o {output.n} 2>> {log}
        """

def get_raw(subsample, p, n):
    check_val("subsample", subsample, bool)
    check_val("p[seqkit]", p, float)
    check_val("n[seqkit]", n, int)
    if subsample == True:
        return rules.subsample.output.n
    else:
        return rules.collect_fastq.output

# trim primers 
rule trim_primers:
    input: get_raw(config["subsample"], config["seqkit"]["p"], config["seqkit"]["n"])
    output: temp(OUTPUT_DIR + "/qc/primers_trimmed/{barcode}.fastq")
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5_patterns,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = config["cutadapt"]["minimum-length"],
    log: OUTPUT_DIR + "/logs/qc/{barcode}/trim_primers.log"
    benchmark: OUTPUT_DIR + "/benchmarks/qc/{barcode}/trim_primers.txt"
    threads: config["threads"]["large"]
    shell:
        """
        cutadapt \
            -j {threads} \
            -e {params.e} -O {params.O} -m {params.m}\
            --discard-untrimmed \
            {params.f} \
            -o {output} \
            {input} \
            > {log} 2>&1
        """

# quality filter
rule q_filter:
    input:  rules.trim_primers.output,
    output: OUTPUT_DIR + "/qc/qfilt/{barcode}.fastq"
    conda: "../envs/seqkit.yaml"
    params:
        Q = config["seqkit"]["min-qual"],
        m = config["seqkit"]["min-len"],
        M = config["seqkit"]["max-len"],
    log: OUTPUT_DIR + "/logs/qc/{barcode}/q_filter.log"
    benchmark: OUTPUT_DIR + "/benchmarks/qc/{barcode}/q_filter.txt"
    threads: config["threads"]["normal"]
    shell: "seqkit seq -j {threads} -Q {params.Q} -m {params.m} -M {params.M} -i {input} > {output} 2> {log}"

checkpoint exclude_empty_fqs:
    input: lambda wc: expand(OUTPUT_DIR + "/qc/qfilt/{barcode}.fastq", barcode=get_demultiplexed(wc))
    output: directory(OUTPUT_DIR + "/qc/qfilt/empty")
    run:
        import shutil
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            if os.stat(i).st_size == 0:
                shutil.move(i, output[0])

def get_qced(wildcards):
    barcodes = get_demultiplexed(wildcards)
    barcodes_empty = glob_wildcards(checkpoints.exclude_empty_fqs.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}.fastq").barcode
    barcodes_empty = sorted(set(barcodes_empty))
    barcodes = [b for b in barcodes if b not in barcodes_empty]
    return barcodes

#  pooling fqs for sensitivity 
rule combine_fastq:
    input: lambda wc: expand(OUTPUT_DIR + "/qc/qfilt/{barcode}.fastq", barcode=get_qced(wc))
    output: temp(OUTPUT_DIR + "/qc/qfilt/pooled.fastq")
    shell:
        "cat {input} > {output}"

def get_filt(wildcards, pooling = True):
    barcodes = get_qced(wildcards) 
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes.append("pooled")
    return expand(OUTPUT_DIR + "/qc/qfilt/{barcode}.fastq", barcode=barcodes)
