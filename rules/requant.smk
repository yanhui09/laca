# check the existence of the requant directory
def requant_check(path, chimera_check = True):
    requant_path = path + '/f2requant'
    if not os.path.exists(requant_path):
        raise ValueError("\n  f2requant directory not found.\n\tMake sure {} is loaded.\n".format(requant_path))
    else:
        fs = [f for f in os.listdir(requant_path) if os.path.isfile(os.path.join(requant_path, f))]
        # file with suffix .fasta, fa, or fna
        fas = [f for f in fs if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.fna')]
        # file with suffix .fastq. .fq
        fqs = [f for f in fs if f.endswith('.fastq')]
        if len(fas) == 0 or len(fqs) == 0:
            raise ValueError("\n  The required fasta (.fasta|.fa|.fna) or fastq (.fastq|.fq) files \n  are not found in {}.\n".format(requant_path))
        else:
            return chimeraF(chimera_check)

checkpoint requant_dir:
    input: requant_check(OUTPUT_DIR, config["chimeraF"])
    output: touch(directory(OUTPUT_DIR + "/f2requant"))

def col_info_requant(wildcards, type):
    if type == "fq":
        f = glob_wildcards(checkpoints.requant_dir.get(**wildcards).output[0]
        + "/{f, (.*).(fastq|fq)}").f
    elif type == "fa":
        f = glob_wildcards(checkpoints.requant_dir.get(**wildcards).output[0]
        + "/{f, (.*).(fasta|fna|fa)}").f
    else:
        raise ValueError("\n  Only fa, fq as valid types.\n")
    f = sorted(set(f))
    return expand(OUTPUT_DIR + "/f2requant/{f}", f=f)

rule concatenate_ref:
    input: lambda wildcards: col_info_requant(wildcards, "fa"),
    output: temp(OUTPUT_DIR + "/requant/combined.fna")
    log: OUTPUT_DIR + "/logs/requant/concatenate_ref.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/concatenate_ref.txt"
    shell:
        "cat {input} > {output} 2> {log}"

# dereplicate sequences with mmseqs
rule derep_denoised_seqs_re:
    input: 
        first = rules.concatenate_ref.output,
    output: 
        rep = temp(OUTPUT_DIR + "/requant/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp(OUTPUT_DIR + "/requant/mmseqs_all_seqs.fasta"),
        tsv = temp(OUTPUT_DIR + "/requant/mmseqs_cluster.tsv"),
        tmp = temp(directory(OUTPUT_DIR + "/tmp")),
    message: "Dereplicate denoised sequences"
    params:
        prefix = OUTPUT_DIR + "/requant/mmseqs",
        mid = config["mmseqs"]["min-seq-id"],
        c = config["mmseqs"]["c"],
    conda: "../envs/mmseqs2.yaml"
    log: OUTPUT_DIR + "/logs/requant/derep_denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/derep_denoised_seqs.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input.first} {params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.mid} -c {params.c} > {log} 2>&1"

# rm duplicates of reverse complements
rule rmdup_revcom_re:
    input: rules.derep_denoised_seqs_re.output.rep
    output: temp(OUTPUT_DIR + "/requant/mmseqs_rep_seq_rmdup.fasta")
    message: "Remove duplicates of reverse complements"
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/requant/rmdup_revcom.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/rmdup_revcom.txt"
    shell: "seqkit rmdup -s {input} -o {output} -w 0 > {log} 2>&1"

# keep fasta header unique
rule rename_fasta_header_re:
    input: rules.rmdup_revcom_re.output
    output: OUTPUT_DIR + "/rep_seqs_requant.fasta"
    log: OUTPUT_DIR + "/logs/requant/rename_fasta.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/rename_fasta.txt"
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
rule index_re:
    input: rules.rename_fasta_header_re.output,
    output: temp(OUTPUT_DIR + "/rep_seqs_requant.mmi")
    message: "Index denoised sequences [Generate abundance matrix]"
    params:
        index_size = config["minimap"]["index_size"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/requant/index.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/index.txt"
    shell: "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule dict_re:
    input: rules.rename_fasta_header_re.output,
    output: temp(OUTPUT_DIR + "/rep_seqs_requant.dict")
    message: "Dict denoised sequences [Generate abundance matrix]"
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/requant/dict.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/dict.txt"
    shell: "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap_rep_seqs_re:
    input:
        fq = OUTPUT_DIR + "/f2requant/{barcode}.fastq",
        mmi = rules.index_re.output,
        dict = rules.dict_re.output,
    output: temp(OUTPUT_DIR + "/requant/mapped/{barcode}.bam")
    message: "Re-map {wildcards.barcode}.fastq files [Generate abundance matrix]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/requant/minimap/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/minimap/{barcode}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        minimap2 -t {threads} -ax {params.x} --secondary=no {input.mmi} {input.fq} 2> {log} | \
        grep -v "^@" | cat {input.dict} - | \
        samtools view -F 3584 -b - > {output} 2>> {log}
        """

rule sort_re:
    input: rules.minimap_rep_seqs_re.output
    output: temp(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam")
    params:
        prefix = OUTPUT_DIR + "/requant/mapped/tmp.{barcode}",
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/requant/sort/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/sort/{barcode}.txt"
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m {params.m} -o {output} 2>{log}"

rule samtools_index_re:        
    input: rules.sort_re.output
    output: temp(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam.bai")
    params:
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/requant/index/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/index/{barcode}.txt"
    shell:
        "samtools index -m {params.m} -@ 1 {input} {output} 2>{log}"

def get_qout_re(wildcards, type_o):
    barcodes = col_info_requant(wildcards, "fq")
    # extract basename & strip extension
    barcodes = [os.path.basename(barcode).split(".")[0] for barcode in barcodes]
    if type_o == "bam":
        output = expand(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand(OUTPUT_DIR + "/requant/mapped/{barcode}.count", barcode=barcodes)
    else:
        raise ValueError("type_o must be 'bam', 'bai', or 'count'")
    return output

# biom format header
rule rowname_kOTU_re:
    input:
        bam = lambda wildcards: get_qout_re(wildcards, "bam"),
        bai = lambda wildcards: get_qout_re(wildcards, "bai"),
    output: temp(OUTPUT_DIR + "/requant/rowname_seqs")
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/requant/rowname_kOTU.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/rowname_kOTU.txt"
    shell:
        """
        echo '#OTU ID' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule seqs_count_re:
    input:
        bam = OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam",
        bai = OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam.bai"
    output: temp(OUTPUT_DIR + "/requant/mapped/{barcode}.count")
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/requant/seqs_count/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/seqs_count/{barcode}.txt"
    shell:
        """
        echo '{wildcards.barcode}' > {output}
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 >> {output}
        """

rule count_matrix_re:
    input:
        rowname_seqs = rules.rowname_kOTU_re.output,
        seqs_count = lambda wildcards: get_qout_re(wildcards, "count"),
    output: OUTPUT_DIR + "/count_matrix_requant.tsv"
    shell:
        "paste {input.rowname_seqs} {input.seqs_count} > {output}"

