# primer info
fprimers = config["fprimer"]
rprimers = config["rprimer"]

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
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

rule subsample:
    input: rules.collect_fastq.output
    output:
        p = temp("qc/subsampled/{barcode}_p.fastq"),
        n = temp("qc/subsampled/{barcode}.fastq"),
    params:
        n = config["seqkit"]["n"],
    log: "logs/qc/subsample/{barcode}.log"
    benchmark: "benchmarks/qc/subsample/{barcode}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        nlines=$(cat {input} | wc -l)
        nreads=$((nlines / 4))
        p=$(echo "scale=1; {params.n} / $nreads + 0.1" | bc)
        if (( $(echo "$p > 1" | bc -l) )); then
            p=1
        fi
        seqkit sample -p $p -j {threads} {input} -o {output.p} -w0 -s123 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -o {output.n} -w0 2>> {log}
        """

def get_raw(subsample = config["subsample"], n = config["seqkit"]["n"]):
    check_val("subsample", subsample, bool)
    check_val("n[seqkit]", n, int)
    if subsample is True:
        return rules.subsample.output.n
    else:
        return rules.collect_fastq.output

# check primer-pattern, process two strands independently
rule check_primers:
    input: get_raw()
    output: 
        passed = temp("qc/primers_passed/{barcode}F.fastq"),
        unpassed = temp("qc/primers_unpassed/{barcode}F.fastq"),
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5_pattern1,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = config["seqkit"]["min_len"],
        M = config["seqkit"]["max_len"],
        action = config["cutadapt"]["action"],
    log: "logs/qc/check_primersF/{barcode}.log"
    benchmark: "benchmarks/qc/check_primersF/{barcode}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        cutadapt \
        --action={params.action} \
        -j {threads} \
        -e {params.e} -O {params.O} -m {params.m}  -M {params.M} \
        {params.f} \
        --untrimmed-output {output.unpassed} \
        -o {output.passed} \
        {input} \
        > {log} 2>&1
        """

use rule check_primers as check_primersR with:
    input: 
        rules.check_primers.output.unpassed
    output:
        passed = temp("qc/primers_passed/{barcode}R.fastq"),
        unpassed = temp("qc/primers_unpassed/{barcode}.fastq"),
    params:
        f = f5_pattern2,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = config["seqkit"]["min_len"],
        M = config["seqkit"]["max_len"],
        action = config["cutadapt"]["action"],
    log: 
        "logs/qc/check_primersR/{barcode}.log"
    benchmark: 
        "benchmarks/qc/check_primersR/{barcode}.txt"

# reverse complement for reverse strand; combine two strands
rule revcomp_fq_combine:
    input: 
        primerF = rules.check_primers.output.passed,
        primerR = rules.check_primersR.output.passed,
    output: 
        revcompR = temp("qc/primers_passed/{barcode}R_revcomp.fastq"),
        combined = temp("qc/primers_passed/{barcode}.fastq"),
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: 
        """
        seqkit seq -j {threads} -r -p -g -t dna {input.primerR} > {output.revcompR} --quiet
        cat {input.primerF} {output.revcompR} > {output.combined}
        """

# option to trim or not
def get_primer_check(primer_check = config["primer_check"], subsample = config["subsample"], n = config["seqkit"]["n"]):
    check_val("primer_check", primer_check, bool)
    if primer_check is True:
        return rules.revcomp_fq_combine.output.combined
    else:
        return get_raw(subsample, n)

# filter chimeric reads
rule minimap2ava_yacrd:
    input: get_primer_check()
    output: temp("qc/yacrd/{barcode}.paf")
    conda: "../envs/yacrd.yaml"
    params:
        x = config["minimap2"]["x_ava"],
        g = config["yacrd"]["minimap2"]["g"],
        f = config["minimap2"]["f"],
    log: "logs/qc/yacrd/{barcode}_ava.log"
    benchmark: "benchmarks/qc/yacrd/{barcode}_ava.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["simple"],
    shell: "minimap2 -x {params.x} -g {params.g} -f {params.f} -t {threads} {input} {input} > {output} 2> {log}"

rule yacrd:
    input: 
        fq = get_primer_check(),
        ava = rules.minimap2ava_yacrd.output
    output: temp("qc/yacrd/{barcode}.fastq")
    conda: "../envs/yacrd.yaml"
    params:
        c = config["yacrd"]["c"],
        n = config["yacrd"]["n"],
    log: "logs/qc/yacrd/{barcode}_filter.log"
    benchmark: "benchmarks/qc/yacrd/{barcode}_filter.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["simple"],
    shell: "yacrd -i {input.ava} -o {log} -c {params.c} -n {params.n} -t {threads} filter -i {input.fq} -o {output} 2>> {log}"

def get_chimera_free(chimera_filt= config["chimera_filt"], primer_check = config["primer_check"], subsample = config["subsample"], n = config["seqkit"]["n"]):
    check_val("chimera_filt", chimera_filt, bool)
    if chimera_filt is True:
        return rules.yacrd.output
    else:
        return get_primer_check(primer_check, subsample, n)

rule q_filter:
    input: get_chimera_free()
    output: "qc/qfilt/{barcode}.fastq"
    params:
        Q = config["seqkit"]["min_qual"],
        m = config["seqkit"]["min_len"],
        M = config["seqkit"]["max_len"],
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: "seqkit seq -j {threads} -Q {params.Q} -m {params.m} -M {params.M} -i -g {input} > {output} --quiet"

checkpoint exclude_empty_fqs:
    input: lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc))
    output: touch(".qc_DONE")

def get_qced_barcodes(wildcards):
    barcodes = get_demux_barcodes(wildcards)
    checkflag = checkpoints.exclude_empty_fqs.get(**wildcards).output[0]
    for i in barcodes:
        if os.stat("qc/qfilt/" + i + ".fastq").st_size == 0:
            barcodes.remove(i)
    return barcodes

localrules: combine_fastq
#  sample pooling to increase sensitivity 
rule combine_fastq:
    input: lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_qced_barcodes(wc))
    output: "qc/qfilt/pooled.fastq"
    shell: "cat {input} > {output}"

def get_filt(wildcards, pool = config["pool"]):
    barcodes = get_qced_barcodes(wildcards) 
    check_val("pool", pool, bool)
    if pool is True:
        barcodes.append("pooled")
    return expand("qc/qfilt/{barcode}.fastq", barcode=barcodes)