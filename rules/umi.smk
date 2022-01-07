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

# indepent qfilt (qc.smk trims the primer together with linker and umi)
use rule q_filter as qfilter_umi with:
    input:
        get_raw(config["subsample"], config["seqkit"]["p"], config["seqkit"]["n"])
    output:
        OUTPUT_DIR + "/umi/{barcode}/qfilt.fastq"
    params:
        Q = config["umi"]["seqkit"]["min-qual"],
        m = config["umi"]["seqkit"]["min-len"],
        M = config["umi"]["seqkit"]["max-len"],
    log:
        OUTPUT_DIR + "/logs/umi/{barcode}/qfilter.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/umi/{barcode}/qfilter.txt"
    threads:
        config["threads"]["large"]

# trim umi region
rule umi_loc:
    input: rules.qfilter_umi.output
    output:
        start=OUTPUT_DIR + "/umi/{barcode}/start.fastq",
        end=OUTPUT_DIR + "/umi/{barcode}/end.fastq",
    conda: "../envs/seqkit.yaml"
    params:
        umi_loc=config["umi"]["loc"]
    log: OUTPUT_DIR + "/logs/umi/{barcode}/umi_loc.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/umi_loc.txt"
    shell:
        """
        seqkit subseq -r 1:{params.umi_loc} {input} 2> {log} 1> {output.start}
        seqkit subseq -r -{params.umi_loc}:-1 {input} 2>> {log} 1> {output.end}
        """

# extract UMI sequences
# linked primers to ensure a "trusted" umi, and take a more trusted ref with adequate size in cluster
# This is different from the original design due to different UMI settings
# https://github.com/SorenKarst/longread_umi/blob/00302fd34cdf7a5b8722965f3f6c581acbafd70c/scripts/umi_binning.sh#L193
rule extract_umi:
    input:
        start=rules.umi_loc.output.start,
        end=rules.umi_loc.output.end,
    output:
        umi1=OUTPUT_DIR + "/umi/{barcode}/umi1.fastq",
        umi2=OUTPUT_DIR + "/umi/{barcode}/umi2.fastq",
    conda: "../envs/cutadapt.yaml"
    params:
        f=f_pattern,
        r=r_pattern,
        max_err=config["umi"]["max_err"],
        min_overlap=config["umi"]["min_overlap"],
        min_len=config["umi"]["len"] - config["umi"]["base_flex"],
        max_len=config["umi"]["len"] + config["umi"]["base_flex"],
    log: OUTPUT_DIR + "/logs/umi/{barcode}/extract_umi.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/extract_umi.txt"
    threads: config["threads"]["normal"]
    shell:
        "cutadapt "
        "-j {threads} -e {params.max_err} -O {params.min_overlap} "
        "-m {params.min_len} -M {params.max_len} "
        "--discard-untrimmed "
        "{params.f} "
        "{params.r} "
        "-o {output.umi1} -p {output.umi2} "
        "{input.start} {input.end} > {log} 2>&1"

# combine UMI sequences
rule concat_umi:
    input:
        umi1=rules.extract_umi.output.umi1,
        umi2=rules.extract_umi.output.umi2,
    output: OUTPUT_DIR + "/umi/{barcode}/umi.fasta"
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/umi/{barcode}/concat_umi.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/concat_umi.txt"
    shell: "seqkit concat {input.umi1} {input.umi2} 2> {log} | seqkit fq2fa -o {output} 2>> {log}"

# cluster UMIs
rule cluster_umi:
    input: rules.concat_umi.output
    output:
        centroid=OUTPUT_DIR + "/umi/{barcode}/centroid.fasta",
        consensus=OUTPUT_DIR + "/umi/{barcode}/consensus.fasta",
        uc=OUTPUT_DIR + "/umi/{barcode}/uc.txt"  
    log: OUTPUT_DIR + "/logs/umi/{barcode}/cluster_umi.log"
    threads: config["threads"]["normal"]
    conda: "../envs/vsearch.yaml"
    params:
        cl_identity = config["umi"]["cl_identity"],
        min_len = 2*(config["umi"]["len"] - config["umi"]["base_flex"]),
        max_len = 2*(config["umi"]["len"] + config["umi"]["base_flex"]),
    shell:
        "vsearch "
        "--cluster_fast {input} --clusterout_sort -id {params.cl_identity} "
        "--clusterout_id --sizeout -uc {output.uc} "
        "--centroids {output.centroid} --consout {output.consensus} "
        "--minseqlength {params.min_len} --maxseqlength {params.max_len} "
        "--qmask none --threads {threads} "
        "--strand both > {log} 2>&1"

rule rename_umi_centroid:
    input: rules.cluster_umi.output.centroid
    output: OUTPUT_DIR + "/umi/{barcode}/umi_centroid.fasta"
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/umi/{barcode}/rename_umi_centroid.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/rename_umi_centroid.txt"
    shell: "seqkit replace -p '^(.+)$' -r 'umi{{nr}}' {input} -o {output} 2> {log}"

# align the trimmed UMI reads to the potential UMIs, which determine the final UMI clusters

def trim_primers(primers, reverse=False):
    if isinstance(primers, str):
        primers = dict(seqs = primers)
    patterns = {k: '-g ' + v for k, v in primers.items()}
    if reverse:
        patterns = {k: '-G ' + v for k, v in primers.items()}
    return ' '.join(v for v in patterns.values())

fprimer1_trim = trim_primers(flinker)
fprimer2_trim = trim_primers(rlinker)
fprimers_trim = fprimer1_trim + ' ' + fprimer2_trim

rprimer1_trim = trim_primers(fprimersR, reverse=True)
rprimer2_trim = trim_primers(rprimersR, reverse=True)
rprimers_trim = rprimer1_trim + ' ' + rprimer2_trim

rule extract_umip:
    input:
        start=rules.umi_loc.output.start,
        end=rules.umi_loc.output.end,
    output:
        umi1=OUTPUT_DIR + "/umi/{barcode}/umi1p.fastq",
        umi2=OUTPUT_DIR + "/umi/{barcode}/umi2p.fastq",
    params:
        f=fprimers_trim,
        r=rprimers_trim,
        max_err=config["umi"]["max_err"],
        min_overlap=config["umi"]["min_overlap"],
        min_len=config["umi"]["len"] - config["umi"]["base_flex"],
        max_len=config["umi"]["len"] + config["umi"]["base_flex"],
    log:
        OUTPUT_DIR + "/logs/umi/{barcode}/extract_umip.log"
    benchmark:
        OUTPUT_DIR + "/benchmarks/umi/{barcode}/extract_umip.txt"
    threads: config["threads"]["normal"]
    shell:
        "cutadapt "
        "-j {threads} -e {params.max_err} -O {params.min_overlap} "
        "-m {params.min_len} -l {params.max_len} "
        "--discard-untrimmed "
        "{params.f} "
        "{params.r} "
        "-o {output.umi1} -p {output.umi2} "
        "{input.start} {input.end} > {log} 2>&1"

rule concat_umip:
    input:
        umi1=rules.extract_umip.output.umi1,
        umi2=rules.extract_umip.output.umi2,
    output: OUTPUT_DIR + "/umi/{barcode}/umip.fastq"
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/umi/{barcode}/concat_umip.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/concat_umip.txt"
    shell: "seqkit concat {input.umi1} {input.umi2} -o {output} 2> {log}"

# calculate the UMI cluster size through mapping
# use bwa aln to map the potential UMI reads to the UMI ref, considering the limited length and sensitivity
rule bwa_index:
    input:
        ref = rules.rename_umi_centroid.output,
    output:
        multiext(OUTPUT_DIR + "/umi/{barcode}/umi_centroid.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    conda: "../envs/bwa.yaml"
    log: OUTPUT_DIR + "/logs/umi/{barcode}/bwa_index.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/bwa_index.txt"
    shell: "bwa index {input.ref} 2> {log}"

rule bwa_aln:
    input:
        rules.bwa_index.output,
        umip = rules.concat_umip.output,
        ref = rules.rename_umi_centroid.output,
    output: OUTPUT_DIR + "/umi/{barcode}/umip.sai",
    conda: "../envs/bwa.yaml"
    log: OUTPUT_DIR + "/logs/umi/{barcode}/bwa_aln.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/bwa_aln.txt"
    threads: config["threads"]["large"]
    shell:
        "bwa aln -n 6 -t {threads} -N {input.ref} {input.umip} > {output} 2> {log}"

rule bwa_samse:
    input:
        umip = rules.concat_umip.output,
        ref = rules.rename_umi_centroid.output,
        sai = rules.bwa_aln.output,
    output: OUTPUT_DIR + "/umi/{barcode}/umip.sam",
    conda: "../envs/bwa.yaml"
    log: OUTPUT_DIR + "/logs/umi/{barcode}/bwa_samse.log"
    benchmark: OUTPUT_DIR + "/benchmarks/umi/{barcode}/bwa_samse.txt"
    shell:
        "bwa samse -n 10000000 {input.ref} {input.sai} {input.umip} 2> {log} | "
        "samtools view -F 4 - > {output} 2>> {log}"