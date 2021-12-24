# check classifier choice
check_list_ele("cluster", config["cluster"], ["clustCon", "isONclustCon", "isONcorCon"])

# dereplicate sequences with mmseqs
rule derep_denoised_seqs:
    input: 
        [OUTPUT_DIR + "/" + str(x) + ".fna" for x in config["cluster"]][1:],
        first = OUTPUT_DIR + "/" + str(config["cluster"][0]) + ".fna",
    output: 
        rep = temp(OUTPUT_DIR + "/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp(OUTPUT_DIR + "/mmseqs_all_seqs.fasta"),
        tsv = temp(OUTPUT_DIR + "/mmseqs_cluster.tsv"),
        tmp = temp(directory(OUTPUT_DIR + "/tmp")),
    message: "Dereplicate denoised sequences"
    params:
        prefix = OUTPUT_DIR + "/mmseqs",
        mid = config["mmseqs"]["min-seq-id"],
        c = config["mmseqs"]["c"],
    conda: "../envs/mmseqs2.yaml"
    log: OUTPUT_DIR + "/logs/derep_denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/derep_denoised_seqs.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input.first} {params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.mid} -c {params.c} > {log} 2>&1"

# keep fasta header unique
rule rename_fasta_header:
    input: rules.derep_denoised_seqs.output.rep
    output: OUTPUT_DIR + "/rep_seqs.fasta"
    log: OUTPUT_DIR + "/logs/rename_fasta.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rename_fasta.benchmark"
    run:
        with open(output[0], "w") as out:
            with open (input[0], "r") as inp:
                i = 1
                for line in inp:
                    if line.startswith(">"):
                        line = ">kOTU_" + str(i) + "\n"
                        i += 1 
                    out.write(line)

# create abundance matrix with minimap
rule index:
    input: rules.rename_fasta_header.output,
    output: temp(OUTPUT_DIR + "/mmseqs_rep_seq.mmi")
    message: "Index denoised sequences [Generate abundance matrix]"
    params:
        index_size = config["minimap"]["index_size"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/index.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/index.txt"
    shell: "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule dict:
    input: rules.rename_fasta_header.output,
    output: temp(OUTPUT_DIR + "/mmseqs_rep_seq.dict")
    message: "Dict denoised sequences [Generate abundance matrix]"
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/dict.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/dict.txt"
    shell: "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap_rep_seqs:
    input:
        fq = rules.q_filter.output,
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
        minimap2 -t {threads} -ax {params.x} --secondary=no {input.mmi} {input.fq} | \
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

def get_qout(wildcards, type_o):
    barcodes = get_demultiplexed(wildcards)
    if type_o == "bam":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand(OUTPUT_DIR + "/mapped/{barcode}.count", barcode=barcodes)
    else:
        raise ValueError("type_o must be 'bam', 'bai', or 'count'")
    return output

# biom format header
rule rowname_kOTU:
    input:
        bam = lambda wildcards: get_qout(wildcards, "bam"),
        bai = lambda wildcards: get_qout(wildcards, "bai"),
    output: temp(OUTPUT_DIR + "/rowname_seqs")
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/rowname_kOTU.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/rowname_kOTU.txt"
    shell:
        """
        echo '#OTU ID' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule seqs_count:
    input:
        bam = OUTPUT_DIR + "/mapped/{barcode}.sorted.bam",
        bai = OUTPUT_DIR + "/mapped/{barcode}.sorted.bam.bai"
    output: temp(OUTPUT_DIR + "/mapped/{barcode}.count")
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/rep_seqs/seqs_count/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/rep_seqs/seqs_count/{barcode}.txt"
    shell:
        """
        echo '{wildcards.barcode}' > {output}
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 >> {output}
        """

rule count_matrix:
    input:
        rowname_seqs = rules.rowname_kOTU.output,
        seqs_count = lambda wildcards: get_qout(wildcards, "count"),
    output: OUTPUT_DIR + "/count_matrix.tsv"
    shell:
        "paste {input.rowname_seqs} {input.seqs_count} > {output}"