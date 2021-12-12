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

# pooling mode will be defined later 
rule combine_fastq:
    input: get_demultiplexed
    output: temp(OUTPUT_DIR + "/raw_fq/pooled.fastq")
    shell:
        "cat {input} > {output}"

# trim primers 
rule trim_primers:
    input: OUTPUT_DIR + "/raw_fq/{barcode}.fastq"
    output: OUTPUT_DIR + "/umap/{barcode}/trimmed.fastq",
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5_patterns,
        max_err = config["max_err"],
        min_overlap = config["min_overlap"],
    log: OUTPUT_DIR + "/logs/umap/{barcode}/trim_primers.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umap/{barcode}/trim_primers.txt"
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
rule nanofilt_umap:
    input:  rules.trim_primers.output,
    output: OUTPUT_DIR + "/umap/{barcode}/filt.fastq"
    conda: "../envs/nanofilt.yaml"
    params:
        q = config["NanoFilt"]["q"],
        l = config["NanoFilt"]["l"],
        maxlength = config["NanoFilt"]["maxlength"],
    log: OUTPUT_DIR + "/logs/umap/{barcode}/nanofilt.log"
    shell: 
        "cat {input} | NanoFilt -q {params.q} -l {params.l} --maxlength {params.maxlength} 2> {log} 1> {output}"

def get_filt(wildcards, pooling = True):
    barcodes = glob_wildcards(checkpoints.guppy_demultiplex.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    barcodes = list(set(barcodes)) 
    if pooling:
        fqs = expand(OUTPUT_DIR + "/umap/{barcode}/trimmed.fastq", barcode=barcodes)
        fqs.append(OUTPUT_DIR + "/umap/pooled/filt.fastq")
    elif pooling is False:
        fqs = expand(OUTPUT_DIR + "/umap/{barcode}/filt.fastq", barcode=barcodes)
    else:
        raise ValueError('Pooling only allows bool type [True/False].\n{} is used in the config file'.format(x))
    return fqs
