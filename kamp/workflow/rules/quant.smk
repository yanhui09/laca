# check classifier choice
check_list_ele("cluster", config["cluster"], ["kmerCon", "clustCon", "isONclustCon", "isONcorCon", "umiCon"])

# dereplicate sequences with mmseqs
rule derep_denoised_seqs:
    input: 
        [OUTPUT_DIR + "/" + str(x) + ".fna" for x in config["cluster"]][1:],
        first = OUTPUT_DIR + "/" + str(config["cluster"][0]) + ".fna",
    output: 
        rep = temp(OUTPUT_DIR + "/quant/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp(OUTPUT_DIR + "/quant/mmseqs_all_seqs.fasta"),
        tsv = temp(OUTPUT_DIR + "/quant/mmseqs_cluster.tsv"),
        tmp = temp(directory(OUTPUT_DIR + "/tmp")),
    message: "Dereplicate denoised sequences"
    params:
        prefix = OUTPUT_DIR + "/quant/mmseqs",
        mid = config["mmseqs"]["min-seq-id"],
        c = config["mmseqs"]["c"],
    conda: "../envs/mmseqs2.yaml"
    log: OUTPUT_DIR + "/logs/quant/derep_denoised_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/derep_denoised_seqs.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input.first} {params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.mid} -c {params.c} > {log} 2>&1"

# rm duplicates of reverse complements
rule rmdup_revcom:
    input: rules.derep_denoised_seqs.output.rep
    output: temp(OUTPUT_DIR + "/quant/mmseqs_rep_seq_rmdup.fasta")
    message: "Remove duplicates of reverse complements"
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/quant/rmdup_revcom.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/rmdup_revcom.txt"
    shell: "seqkit rmdup -s {input} -o {output} -w 0 > {log} 2>&1"

# keep fasta header unique
rule rename_fasta_header:
    input: rules.rmdup_revcom.output
    output: OUTPUT_DIR + "/rep_seqs.fasta"
    log: OUTPUT_DIR + "/logs/quant/rename_fasta.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/rename_fasta.benchmark"
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
    output: temp(OUTPUT_DIR + "/rep_seqs.mmi")
    message: "Index denoised sequences [Generate abundance matrix]"
    params:
        index_size = config["minimap"]["index_size"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/quant/index.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/index.txt"
    shell: "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule dict:
    input: rules.rename_fasta_header.output,
    output: temp(OUTPUT_DIR + "/rep_seqs.dict")
    message: "Dict denoised sequences [Generate abundance matrix]"
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/quant/dict.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/dict.txt"
    shell: "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap_rep_seqs:
    input:
        fq = rules.q_filter.output,
        mmi = rules.index.output,
        dict = rules.dict.output,
    output: temp(OUTPUT_DIR + "/quant/mapped/{barcode}.bam")
    message: "Re-map {wildcards.barcode}.fastq files [Generate abundance matrix]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/quant/minimap/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/minimap/{barcode}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        minimap2 -t {threads} -ax {params.x} --secondary=no {input.mmi} {input.fq} 2> {log} | \
        grep -v "^@" | cat {input.dict} - | \
        samtools view -F 3584 -b - > {output} 2>> {log}
        """

rule sort:
    input: rules.minimap_rep_seqs.output
    output: temp(OUTPUT_DIR + "/quant/mapped/{barcode}.sorted.bam")
    params:
        prefix = OUTPUT_DIR + "/quant/mapped/tmp.{barcode}",
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/quant/sort/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/sort/{barcode}.txt"
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m {params.m} -o {output} 2>{log}"

rule samtools_index:        
    input: rules.sort.output
    output: temp(OUTPUT_DIR + "/quant/mapped/{barcode}.sorted.bam.bai")
    params:
        m = config["samtools"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/quant/index/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/index/{barcode}.txt"
    shell:
        "samtools index -m {params.m} -@ 1 {input} {output} 2>{log}"

def get_qout(wildcards, type_o):
    barcodes = get_qced(wildcards)
    if type_o == "bam":
        output = expand(OUTPUT_DIR + "/quant/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand(OUTPUT_DIR + "/quant/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand(OUTPUT_DIR + "/quant/mapped/{barcode}.count", barcode=barcodes)
    else:
        raise ValueError("type_o must be 'bam', 'bai', or 'count'")
    return output

# biom format header
rule rowname_kOTU:
    input:
        bam = lambda wildcards: get_qout(wildcards, "bam"),
        bai = lambda wildcards: get_qout(wildcards, "bai"),
    output: temp(OUTPUT_DIR + "/quant/rowname_seqs")
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/quant/rowname_kOTU.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/rowname_kOTU.txt"
    shell:
        """
        echo '#OTU ID' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule seqs_count:
    input:
        bam = OUTPUT_DIR + "/quant/mapped/{barcode}.sorted.bam",
        bai = OUTPUT_DIR + "/quant/mapped/{barcode}.sorted.bam.bai"
    output: temp(OUTPUT_DIR + "/quant/mapped/{barcode}.count")
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/quant/seqs_count/{barcode}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/quant/seqs_count/{barcode}.txt"
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

rule q2_repseqs_import:
    input: rules.rename_fasta_header.output
    output: temp(OUTPUT_DIR + "/rep_seqs.qza")
    conda: "../envs/q2_plugins.yaml"
    log: OUTPUT_DIR + "/logs/chimeraF/q2_repseqs_import.log"
    benchmark: OUTPUT_DIR + "/benchmarks/chimeraF/q2_repseqs_import.txt"
    shell:
        """
        qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path {input} \
        --output-path {output} \
        > {log} 2>&1
        """

rule q2_ftable_import:
    input: rules.count_matrix.output
    output: 
        biom = temp(OUTPUT_DIR + "/table.biom"),
        qza = temp(OUTPUT_DIR + "/table.qza")
    conda: "../envs/q2_plugins.yaml"
    log: OUTPUT_DIR + "/logs/chimeraF/q2_ftable_import.log"
    benchmark: OUTPUT_DIR + "/benchmarks/chimeraF/q2_ftable_import.txt"
    shell:
        """
        biom convert -i {input} -o {output.biom} --table-type="OTU table" --to-hdf5 > {log} 2>&1
        qiime tools import \
        --input-path {output.biom} \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --output-path {output.qza} \
        >> {log} 2>&1
        """

rule q2_uchime_denovo:
    input:
        ftable = rules.q2_ftable_import.output.qza,
        repseqs = rules.q2_repseqs_import.output,
    output: 
        nonchimeras = OUTPUT_DIR + "/chimeraF/nonchimeras.qza",
        chimeras = OUTPUT_DIR + "/chimeraF/chimeras.qza",
        stats = OUTPUT_DIR + "/chimeraF/stats.qza",
    conda: "../envs/q2_plugins.yaml"
    log: OUTPUT_DIR + "/logs/chimeraF/uchime_denovo.log"
    benchmark: OUTPUT_DIR + "/benchmarks/chimeraF/uchime_denovo.txt"
    shell:
        """
        qiime vsearch uchime-denovo \
        --i-table {input.ftable} \
        --i-sequences {input.repseqs} \
        --o-nonchimeras {output.nonchimeras} \
        --o-chimeras {output.chimeras} \
        --o-stats {output.stats} \
        > {log} 2>&1
        """

rule q2_filter_features:
    input:
        ftable = rules.q2_ftable_import.output.qza,
        nonchimeras = rules.q2_uchime_denovo.output.nonchimeras,
    output: temp(OUTPUT_DIR + "/chimeraF/table_nonchimeras.qza")
    conda: "../envs/q2_plugins.yaml"
    log: OUTPUT_DIR + "/logs/chimeraF/filter_features.log"
    benchmark: OUTPUT_DIR + "/benchmarks/chimeraF/filter_features.txt"
    shell:
        """
        qiime feature-table filter-features \
        --i-table {input.ftable} \
        --m-metadata-file {input.nonchimeras} \
        --o-filtered-table {output} \
        > {log} 2>&1
        """

rule q2_filter_seqs:
    input:
        repseqs = rules.q2_repseqs_import.output,
        nonchimeras = rules.q2_uchime_denovo.output.nonchimeras,
    output: temp(OUTPUT_DIR + "/chimeraF/rep_seqs_nonchimeras.qza")
    conda: "../envs/q2_plugins.yaml"
    log: OUTPUT_DIR + "/logs/chimeraF/filter_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/chimeraF/filter_seqs.txt"
    shell:
        """
        qiime feature-table filter-seqs \
        --i-data {input.repseqs} \
        --m-metadata-file {input.nonchimeras} \
        --o-filtered-data {output} \
        > {log} 2>&1
        """

rule q2_ftable_export:
    input: rules.q2_filter_features.output
    output:
        biom = temp(OUTPUT_DIR + "/chimeraF/table_nonchimeras.biom"), 
        tsv = OUTPUT_DIR + "/chimeraF/table_nonchimeras.tsv"
    params:
        _dir = OUTPUT_DIR + "/chimeraF"
    conda: "../envs/q2_plugins.yaml"
    log: OUTPUT_DIR + "/logs/chimeraF/ftable_export.log"
    benchmark: OUTPUT_DIR + "/benchmarks/chimeraF/ftable_export.txt"
    shell:
        """
        qiime tools export \
        --input-path {input} \
        --output-path {params._dir} \
        > {log} 2>&1
        mv {params._dir}/feature-table.biom {output.biom}
        biom convert -i {output.biom} -o {output.tsv} --to-tsv >> {log} 2>&1
        """

rule q2_repseqs_export:
    input: rules.q2_filter_seqs.output
    output: OUTPUT_DIR + "/chimeraF/rep_seqs_nonchimeras.fasta"
    conda: "../envs/q2_plugins.yaml"
    params:
        _dir = OUTPUT_DIR + "/chimeraF"
    log: OUTPUT_DIR + "/logs/chimeraF/repseqs_export.log"
    benchmark: OUTPUT_DIR + "/benchmarks/chimeraF/repseqs_export.txt"
    shell:
        """
        qiime tools export \
        --input-path {input} \
        --output-path {params._dir} \
        > {log} 2>&1
        mv {params._dir}/dna-sequences.fasta {output}
        """

def chimeraF(chimera_check = True):
    check_val("chimeraF", chimera_check, bool)
    fs = ["chimeraF/" + x for x in ["table_nonchimeras.tsv", "rep_seqs_nonchimeras.fasta"]]
    if chimera_check == False:
        fs = ["count_matrix.tsv", "rep_seqs.fasta"]
    return expand(OUTPUT_DIR + "/{f}", f = fs)