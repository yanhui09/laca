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
#---------

rule subsample:
    input: rules.collect_fastq.output
    output:
        p = temp("qc/subsampled/{barcode}_p.fastq"),
        n = temp("qc/subsampled/{barcode}.fastq"),
    conda: "../envs/seqkit.yaml"
    params:
        p = config["seqkit"]["p"],
        n = config["seqkit"]["n"],
    log: "logs/subsample/{barcode}.log"
    benchmark: "benchmarks/subsample/{barcode}.txt"
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
    if subsample is True:
        return rules.subsample.output.n
    else:
        return rules.collect_fastq.output

# trim primers 
# process two strands differently
rule trim_primers:
    input: get_raw(config["subsample"], config["seqkit"]["p"], config["seqkit"]["n"])
    output: 
        trimmed = temp("qc/primers_trimmed/{barcode}F.fastq"),
        untrimmed = temp("qc/primers_untrimmed/{barcode}F.fastq"),
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5_pattern1,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = config["cutadapt"]["minimum-length"],
    log: "logs/qc/{barcode}/trim_primersF.log"
    benchmark: "benchmarks/qc/{barcode}/trim_primersF.txt"
    threads: config["threads"]["large"]
    shell:
        """
        cutadapt \
        -j {threads} \
        -e {params.e} -O {params.O} -m {params.m} \
        {params.f} \
        --untrimmed-output {output.untrimmed} \
        -o {output.trimmed} \
        {input} \
        > {log} 2>&1
        """

use rule trim_primers as trim_primersR with:
    input: 
        rules.trim_primers.output.untrimmed
    output:
        trimmed = temp("qc/primers_trimmed/{barcode}R.fastq"),
        untrimmed = "qc/primers_untrimmed/{barcode}.fastq",
    params:
        f = f5_pattern2,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = config["cutadapt"]["minimum-length"],
    log: 
        "logs/qc/{barcode}/trim_primersR.log"
    benchmark: 
        "benchmarks/qc/{barcode}/trim_primersR.txt"

# reverse complement for reverse strand
rule revcomp_fq:
    input: rules.trim_primersR.output.trimmed
    output: temp("qc/primers_trimmed/{barcode}R_revcomp.fastq")
    conda: "../envs/seqkit.yaml"
    log: "logs/qc/{barcode}/revcomp_fq.log"
    benchmark: "benchmarks/qc/{barcode}/revcomp_fq.txt"
    threads: config["threads"]["normal"]
    shell: "seqkit seq -j {threads} -r -p -t dna {input} > {output} 2> {log}"

# trim primers or not
def trim_check(trim, subsample, p, n):
    check_val("trim", trim, bool)
    out = [rules.trim_primers.output.trimmed, rules.revcomp_fq.output]
    if trim is False:
        out = get_raw(subsample, p, n)
    return out

# quality filter
rule q_filter:
    input:
        trim_check(config["trim"], config["subsample"], config["seqkit"]["p"], config["seqkit"]["n"])
    output: temp("qc/qfilt/{barcode}.fastq")
    conda: "../envs/seqkit.yaml"
    params:
        Q = config["seqkit"]["min-qual"],
        m = config["seqkit"]["min-len"],
        M = config["seqkit"]["max-len"],
    log: "logs/qc/{barcode}/q_filter.log"
    benchmark: "benchmarks/qc/{barcode}/q_filter.txt"
    threads: config["threads"]["normal"]
    shell: "cat {input} | seqkit seq -j {threads} -Q {params.Q} -m {params.m} -M {params.M} -i > {output} 2> {log}"

checkpoint exclude_empty_fqs:
    input: lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demultiplexed(wc))
    output: directory("qc/qfilt/empty")
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

#  pool fqs for sensitivity 
rule combine_fastq:
    input: lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_qced(wc))
    output: temp("qc/qfilt/pooled.fastq")
    shell:
        "cat {input} > {output}"

def get_filt(wildcards, pool = True):
    barcodes = get_qced(wildcards) 
    check_val("pool", pool, bool)
    if pool is True:
        barcodes.append("pooled")
    return expand("qc/qfilt/{barcode}.fastq", barcode=barcodes)
