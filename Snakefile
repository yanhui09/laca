#--------------
configfile: "config.yaml"
wildcard_constraints:
        barcode="BRK[0-9][0-9]"
#--------------
INPUT_DIR = config["basecalled_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
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
#---------------

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all:
    input:
        OUTPUT_DIR + "/ONT-L7-GG.txt"
    
checkpoint demultiplex:
    input: INPUT_DIR
    output: directory(OUTPUT_DIR + "/demultiplexed")
    log: OUTPUT_DIR + "/logs/demultiplex.log"
    threads: config["threads"]["large"]
    params:
        guppy=config["guppy"],
        barcode_kits=config["barcode_kits"],
        use_cuda=config["use_cuda"],
    shell:
        """
        {params.guppy}/guppy_barcoder -i {input} -s {output} -t {threads} --barcode_kits {params.barcode_kits} 2>{log}
        if {params.use_cuda}; then
            {params.guppy}/guppy_barcoder -i {input} -s {output} -t {threads} --barcode_kits {params.barcode_kits} -x auto 2>{log}
        fi
        """

rule nanofilt:
    input:  OUTPUT_DIR + "/demultiplexed/{barcode}"
    output: OUTPUT_DIR + "/filt_fq/{barcode}.fastq"
    log: OUTPUT_DIR + "/logs/nanofilt/{barcode}.log"
    conda: "envs/nanofilt.yaml"
    threads: 1
    shell: 
        """
        cat {input}/*.fastq | NanoFilt -q 8 -l 1000 --maxlength 1600 --headcrop 15 --tailcrop 15 2> {log} 1> {output}
        """

# trim umi region
rule umi_loc:
    input: OUTPUT_DIR + "/filt_fq/{barcode}.fastq"
    output:
        start=OUTPUT_DIR + "/umi/{barcode}/start.fastq",
        end=OUTPUT_DIR + "/umi/{barcode}/end.fastq",
    log: OUTPUT_DIR + "/logs/umi/umi_loc/{barcode}.log"
    conda: "envs/seqkit.yaml"
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
    conda: "envs/cutadapt.yaml"
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
            > {log} 2>&1 &
        """

# combine UMI sequences
rule concat_umi:
    input:
        umi1=rules.extract_umi.output.umi1,
        umi2=rules.extract_umi.output.umi2,
    output: OUTPUT_DIR + "/umi/{barcode}/umi.fastq"
    log: OUTPUT_DIR + "/logs/umi/concat_umi/{barcode}.log"
    threads: 1
    conda: "envs/seqkit.yaml"
    shell:
        "seqkit concat {input.umi1} {input.umi2} 2> {log} > {output}"