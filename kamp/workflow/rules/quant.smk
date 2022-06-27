# check classifier choice
check_list_ele("cluster", config["cluster"], ["kmerCon", "clustCon", "isONclustCon", "isONcorCon", "umiCon"])

# dereplicate sequences with mmseqs
rule derep_denoised_seqs:
    input: 
        ["" + str(x) + ".fna" for x in config["cluster"]][1:],
        first = "" + str(config["cluster"][0]) + ".fna",
    output: 
        rep = temp("quant/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp("quant/mmseqs_all_seqs.fasta"),
        tsv = temp("quant/mmseqs_cluster.tsv"),
        tmp = temp(directory("tmp")),
    message: "Dereplicate denoised sequences"
    params:
        prefix = "quant/mmseqs",
        mid = config["mmseqs"]["min-seq-id"],
        c = config["mmseqs"]["c"],
    conda: "../envs/mmseqs2.yaml"
    log: "logs/quant/derep_denoised_seqs.log"
    benchmark: "benchmarks/quant/derep_denoised_seqs.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input.first} {params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.mid} -c {params.c} > {log} 2>&1"

# rm duplicates of reverse complements
rule rmdup_revcom:
    input: rules.derep_denoised_seqs.output.rep
    output: temp("quant/mmseqs_rep_seq_rmdup.fasta")
    message: "Remove duplicates of reverse complements"
    conda: "../envs/seqkit.yaml"
    log: "logs/quant/rmdup_revcom.log"
    benchmark: "benchmarks/quant/rmdup_revcom.txt"
    shell: "seqkit rmdup -s {input} -o {output} -w 0 > {log} 2>&1"

# keep fasta header unique
rule rename_fasta_header:
    input: rules.rmdup_revcom.output
    output: "rep_seqs.fasta"
    log: "logs/quant/rename_fasta.log"
    benchmark: "benchmarks/quant/rename_fasta.benchmark"
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
    output: temp("rep_seqs.mmi")
    message: "Index denoised sequences [Generate abundance matrix]"
    params:
        index_size = config["minimap"]["index_size"],
    conda: "../envs/minimap2.yaml"
    log: "logs/quant/index.log"
    benchmark: "benchmarks/quant/index.txt"
    shell: "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule dict:
    input: rules.rename_fasta_header.output,
    output: temp("rep_seqs.dict")
    message: "Dict denoised sequences [Generate abundance matrix]"
    conda: "../envs/samtools.yaml"
    log: "logs/quant/dict.log"
    benchmark: "benchmarks/quant/dict.txt"
    shell: "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap_rep_seqs:
    input:
        fq = rules.q_filter.output,
        mmi = rules.index.output,
        dict = rules.dict.output,
    output: temp("quant/mapped/{barcode}.bam")
    message: "Re-map {wildcards.barcode}.fastq files [Generate abundance matrix]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/minimap2.yaml"
    log: "logs/quant/minimap/{barcode}.log"
    benchmark: "benchmarks/quant/minimap/{barcode}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        minimap2 -t {threads} -ax {params.x} --secondary=no {input.mmi} {input.fq} 2> {log} | \
        grep -v "^@" | cat {input.dict} - | \
        samtools view -F 3584 -b - > {output} 2>> {log}
        """

rule sort:
    input: rules.minimap_rep_seqs.output
    output: temp("quant/mapped/{barcode}.sorted.bam")
    params:
        prefix = "quant/mapped/tmp.{barcode}",
        m = config["samtools"]["m"],
    conda: "../envs/samtools.yaml"
    log: "logs/quant/sort/{barcode}.log"
    benchmark: "benchmarks/quant/sort/{barcode}.txt"
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m {params.m} -o {output} 2>{log}"

rule samtools_index:        
    input: rules.sort.output
    output: temp("quant/mapped/{barcode}.sorted.bam.bai")
    params:
        m = config["samtools"]["m"],
    conda: "../envs/samtools.yaml"
    log: "logs/quant/index/{barcode}.log"
    benchmark: "benchmarks/quant/index/{barcode}.txt"
    shell:
        "samtools index -m {params.m} -@ 1 {input} {output} 2>{log}"

def get_qout(wildcards, type_o):
    barcodes = get_qced(wildcards)
    if type_o == "bam":
        output = expand("quant/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand("quant/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand("quant/mapped/{barcode}.count", barcode=barcodes)
    else:
        raise ValueError("type_o must be 'bam', 'bai', or 'count'")
    return output

# biom format header
rule rowname_kOTU:
    input:
        bam = lambda wildcards: get_qout(wildcards, "bam"),
        bai = lambda wildcards: get_qout(wildcards, "bai"),
    output: temp("quant/rowname_seqs")
    conda: "../envs/samtools.yaml"
    log: "logs/quant/rowname_kOTU.log"
    benchmark: "benchmarks/quant/rowname_kOTU.txt"
    shell:
        """
        echo '#OTU ID' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule seqs_count:
    input:
        bam = "quant/mapped/{barcode}.sorted.bam",
        bai = "quant/mapped/{barcode}.sorted.bam.bai"
    output: temp("quant/mapped/{barcode}.count")
    conda: "../envs/samtools.yaml"
    log: "logs/quant/seqs_count/{barcode}.log"
    benchmark: "benchmarks/quant/seqs_count/{barcode}.txt"
    shell:
        """
        echo '{wildcards.barcode}' > {output}
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 >> {output}
        """

rule count_matrix:
    input:
        rowname_seqs = rules.rowname_kOTU.output,
        seqs_count = lambda wildcards: get_qout(wildcards, "count"),
    output: "count_matrix.tsv"
    shell:
        "paste {input.rowname_seqs} {input.seqs_count} > {output}"

rule q2_repseqs_import:
    input: rules.rename_fasta_header.output
    output: temp("rep_seqs.qza")
    conda: "../envs/q2_plugins.yaml"
    log: "logs/chimeraF/q2_repseqs_import.log"
    benchmark: "benchmarks/chimeraF/q2_repseqs_import.txt"
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
        biom = temp("table.biom"),
        qza = temp("table.qza")
    conda: "../envs/q2_plugins.yaml"
    log: "logs/chimeraF/q2_ftable_import.log"
    benchmark: "benchmarks/chimeraF/q2_ftable_import.txt"
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
        nonchimeras = "chimeraF/nonchimeras.qza",
        chimeras = "chimeraF/chimeras.qza",
        stats = "chimeraF/stats.qza",
    conda: "../envs/q2_plugins.yaml"
    log: "logs/chimeraF/uchime_denovo.log"
    benchmark: "benchmarks/chimeraF/uchime_denovo.txt"
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
    output: temp("chimeraF/table_nonchimeras.qza")
    conda: "../envs/q2_plugins.yaml"
    log: "logs/chimeraF/filter_features.log"
    benchmark: "benchmarks/chimeraF/filter_features.txt"
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
    output: temp("chimeraF/rep_seqs_nonchimeras.qza")
    conda: "../envs/q2_plugins.yaml"
    log: "logs/chimeraF/filter_seqs.log"
    benchmark: "benchmarks/chimeraF/filter_seqs.txt"
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
        biom = temp("chimeraF/table_nonchimeras.biom"), 
        tsv = "chimeraF/table_nonchimeras.tsv"
    params:
        _dir = "chimeraF"
    conda: "../envs/q2_plugins.yaml"
    log: "logs/chimeraF/ftable_export.log"
    benchmark: "benchmarks/chimeraF/ftable_export.txt"
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
    output: "chimeraF/rep_seqs_nonchimeras.fasta"
    conda: "../envs/q2_plugins.yaml"
    params:
        _dir = "chimeraF"
    log: "logs/chimeraF/repseqs_export.log"
    benchmark: "benchmarks/chimeraF/repseqs_export.txt"
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
    return expand("{f}", f = fs)