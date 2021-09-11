#--------------
configfile: "config.yaml"
wildcard_constraints:
        barcode="BRK[0-9][0-9]"

#--------------
INPUT_DIR = config["basecalled_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
LIVE_BATCH = glob_wildcards(INPUT_DIR + "/{batchid}.fastq").batchid
PROBE_1 = config["probe1_fa"].rstrip("/"),
PROBE_2 = config["probe2_fa"].rstrip("/"),
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
