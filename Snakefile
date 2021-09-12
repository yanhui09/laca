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
        INPUT_DIR
    output: 
        directory(OUTPUT_DIR + "/demultiplexed")
    logs: 
        OUTPUT_DIR + "/logs/demultiplex.log"
    threads: config["threads_demultiplex"]
    params:
        guppy=config["guppy"]
        barcode_kits=config["barcode_kits"],
        use_cuda=config["use_cuda"]
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
    logs: OUTPUT_DIR + "/logs/nanofilt/{barcode}.log"
    conda: "envs/nanofilt.yaml"
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
    logs: OUTPUT_DIR + "/logs/umi_loc/{barcode}.log"
    conda: "envs/seqkit.yaml"
    params:
        umi_loc=config["umi_loc"]
    shell:
        """
        seqkit subseq -r 1:{params.umi_loc} {input} 2> {log} 1> {output.start}
        seqkit subseq -r -{params.umi_loc}:-1 {input} 2>> {log} 1> {output.end}
        """

# extract UMI seqeunces 
rule extract_umi:
    input:
    output:
    logs:
    conda:
        "envs/cutadapt.yaml"
    shell:
        """"

        """"