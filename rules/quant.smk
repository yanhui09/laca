# dereplicate denoised sequences with mmseqs
rule dereplicate_denoised_seqs:
    input: rules.collect_isONclust_consensus.output.fna,
    output: OUTPUT_DIR + "/rep_seqs.fna",
    message: "Dereplicate denoised sequences"
    params:
        tmp = OUTPUT_DIR + "/tmp",
    conda: "../envs/mmseqs2.yaml"
    log: OUTPUT_DIR + "/logs/dereplicate_denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/dereplicate_denoised_seqs.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input} {output} {params.tmp} --threads {threads} --min-seq-id 1 -c 1 > {log} 2>&1"

# create abundance matrix with minimap
rule index:
    input: rules.dereplicate_denoised_seqs.output,
    output: temp(OUTPUT_DIR + "/rep_seqs.mmi")
    message: "Index denoised sequences [Generate abundance matrix]"
    params:
        index_size = config["minimap"]["index_size"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/index.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/index.txt"
    shell: "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule dict:
    input: rules.dereplicate_denoised_seqs.output,
    output: temp(OUTPUT_DIR + "/rep_seqs.dict")
    message: "Dict denoised sequences [Generate abundance matrix]"
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/dict.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/dict.txt"
    shell: "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap_rep_seqs:
    input:
        fq = rules.trim_primers.output,
        mmi = rules.index.output,
        dict = rules.dict.output,
    output: temp(OUTPUT_DIR + "/mapped/{barcode}.bam")
    message: "Re-map trimmed fastq files [Generate abundance matrix]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/minimap/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/minimap/{barcode}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        minimap2 -t {threads} -ax {params.x} {input.mmi} {input.fq} | \
        grep -v "^@" | cat {input.dict} - | \
        samtools view -F 3584 -b - > {output} 2>{log}
        """

rule sort:
    input: rules.minimap_rep_seqs.output
    output: temp(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam")
    params:
        prefix = OUTPUT_DIR + "/mapped/tmp.{barcode}",
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/sort/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/sort/{barcode}.txt"
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m {params.m} -o {output} 2>{log}"

rule samtools_index:        
    input: rules.sort.output
    output: temp(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam.bai")
    params:
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/index/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/index/{barcode}.txt"
    shell:
        "samtools index -m {params.m} -@ 1 {input} {output} 2>{log}"

def get_barcodes(wildcards, type_o):
    barcodes = glob_wildcards(checkpoints.guppy_demultiplex.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    if type_o == "bam":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.count", barcode=barcodes)
    else:
        output = barcodes
    return output

# biom format header
rule header_sample:
    input:
        bai = lambda wildcards: get_barcodes(wildcards, "bai"),
    output: temp(OUTPUT_DIR + "/header_sample")
    run:
        with open(output[0], 'w') as f:
            f.write('#OTU ID\t'+ '\t'.join(SAMPLE) + '\n')

rule rowname_seqs:
    input:
        bam = lambda wildcards: get_barcodes(wildcards, "bam"),
        bai = lambda wildcards: get_barcodes(wildcards, "bai"),
    output: temp(OUTPUT_DIR + "/rowname_seqs")
    conda: "../envs/polish.yaml"
    shell:
        """
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 > {output}
        """
       
rule seqs_count:
    input:
        bam = OUTPUT_DIR + "/mapped/{sample}.sorted.bam",
        bai = OUTPUT_DIR + "/mapped/{sample}.sorted.bam.bai"
    output: temp(OUTPUT_DIR + "/mapped/{sample}.count")
    conda: "../envs/polish.yaml"
    shell:
        """
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 > {output}
        """

rule count_matrix:
    input:
        seqs_count = lambda wildcards: get_barcodes(wildcards, "count"),
        header_sample = rules.header_sample.output,
        rowname_seqs = rules.rowname_seqs.output,
    output: OUTPUT_DIR + "/count_matrix.tsv"
    shell:
        "paste {input.rowname_seqs} {input.seqs_count} | cat {input.header_sample} - > {output}"