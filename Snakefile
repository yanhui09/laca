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
    input:  
        INPUT_DIR + "/{live_batch}.fastq"
    output: 
        directory(OUTPUT_DIR + "/demultiplexed_fq/{live_batch}")
    threads: config["threads_demultiplex"]
    shell:
        """
        mkdir -p {output} && cp {input} {output}
        if nvidia-smi; then
            {config[guppy_dir]}/guppy_barcoder -i {output} -s {output} -t {threads} --barcode_kits SQK16S-GXO -x auto
        else
            {config[guppy_dir]}/guppy_barcoder -i {output} -s {output} -t {threads} --barcode_kits SQK16S-GXO
        fi
        """

rule nanofilt:
    input: 
        OUTPUT_DIR + "/demultiplexed_fq/{live_batch}/{barcode}"
    output: 
        OUTPUT_DIR + "/filt_fq/{live_batch}/{barcode}.fastq"
    conda:
        "envs/nanofilt.yaml"
    shell: 
        """
        cat {input}/*.fastq | NanoFilt -q 8 -l 1000 --maxlength 1600 --headcrop 15 --tailcrop 15 > {output}
        """

# extract UMI sequences (https://github.com/fhlab/UMIC-seq)
rule umi_extract:
    input: 
        rules.nanofilt.output
    output: 
        umi1 = OUTPUT_DIR + "/filt_fq/{barcode}_umi1.fasta",
        umi2 = OUTPUT_DIR + "/filt_fq/{barcode}_umi2.fasta",
    conda:
        "envs/umic-seq.yaml"
    params:
        scripts = "scripts/UMIC-seq.py",
        probe1 = PROBE_1,
        probe2 = PROBE_2,
    shell: 
        """
        python {params.scripts} UMIextract --input {input} --probe {params.probe1} --umi_loc down --umi_len 15 --output {output.umi1}
        python {params.scripts} UMIextract --input {input} --probe {params.probe2} --umi_loc up --umi_len 15 --output {output.umi2}
        """
