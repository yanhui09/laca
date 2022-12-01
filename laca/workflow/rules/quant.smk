# check classifier choice
check_list_ele("cluster", config["cluster"], ["kmerCon", "clustCon", "isONclustCon", "isONcorCon", "umiCon"])

def get_consensus2derep(cls):
   #{cls}/{cls}.fna except umiCon ({cls}/{cls}_trimmed.fna)
    if cls == "umiCon":
        return "{cls}/{cls}_trimmed.fna".format(cls=cls)
    else:
        return "{cls}/{cls}.fna".format(cls=cls)

# dereplicate sequences with MMseqs2
rule derep_denoised_seqs:
    input: 
        [get_consensus2derep(cls=i) for i in config["cluster"]][1:],
        first = get_consensus2derep(cls = config["cluster"][0]),
    output: 
        rep = temp("quant/mmseqs2_rep_seq.fasta"),
        all_by_cluster = temp("quant/mmseqs2_all_seqs.fasta"),
        tsv = temp("quant/mmseqs2_cluster.tsv"),
        tmp = temp(directory("tmp")),
    message: "Dereplicate denoised sequences"
    params:
        prefix = "quant/mmseqs2",
        mid = config["mmseqs2"]["min_id"],
        c = config["mmseqs2"]["c"],
    conda: "../envs/mmseqs2.yaml"
    log: "logs/quant/derep_denoised_seqs.log"
    benchmark: "benchmarks/quant/derep_denoised_seqs.txt"
    threads: config["threads"]["large"]
    shell:
        #https://github.com/soedinglab/MMseqs2/wiki#how-do-parameters-of-cd-hit-relate-to-mmseqs2
        #divergent amplicons for local alignment, mmseqs2 yes, cd-hit no.
        #--min-seq-id 0.99 -c 0.5 (best in benchmark, -c 0.5-0.7 is good)
        "mmseqs easy-cluster {input.first} {params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.mid} -c {params.c} --cluster-reassign > {log} 2>&1"

# rm duplicates of reverse complements
rule rmdup_revcom:
    input: rules.derep_denoised_seqs.output.rep
    output: temp("quant/mmseqs2_rep_seq_rmdup.fasta")
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

# create abundance matrix with minimap2
rule index_repseqs:
    input: rules.rename_fasta_header.output
    output: 
        mmi = temp("rep_seqs.mmi"),
        dict = temp("rep_seqs.dict"),
    message: "Index denoised sequences [Abundance matrix]"
    params:
        index_size = "4G",
    conda: "../envs/minimap2.yaml"
    log: "logs/quant/index_repseqs.log"
    benchmark: "benchmarks/quant/index_repseqs.txt"
    shell: 
        """
        minimap2 -I {params.index_size} -d {output.mmi} {input} 2> {log}
        samtools dict {input} | cut -f1-3 > {output.dict} 2>> {log}
        """

rule minimap2repseqs:
    input:
        fq = rules.q_filter.output,
        mmi = rules.index_repseqs.output.mmi,
        dict = rules.index_repseqs.output.dict,
    output: 
        bam = temp("quant/mapped/{barcode}.bam"),
        sort = temp("quant/mapped/{barcode}.sorted.bam"),
        bai = temp("quant/mapped/{barcode}.sorted.bam.bai"),
        counts = temp("quant/mapped/{barcode}.count"),
    message: "Map reads back, barcode={wildcards.barcode} [Abundance matrix]"
    params:
        x = config["minimap2"]["x_map"],
        prefix = "quant/mapped/tmp.{barcode}",
        m = "3G",
        # https://lh3.github.io/minimap2/minimap2.html#10
        max_de = 0.1,
    conda: "../envs/minimap2.yaml"
    log: "logs/quant/minimap2/{barcode}.log"
    benchmark: "benchmarks/quant/minimap2/{barcode}.txt"
    threads: config["threads"]["normal"]
    shell:
        # de:f: divergence < max_de
        # no supplementary alignments, primary alignments only, exclude unmapped reads
        """
        minimap2 -t {threads} -ax {params.x} {input.mmi} {input.fq} 2> {log} | \
        grep -v "^@" | awk -F '\t|de:f:' '$(NF-1) < {params.max_de}' | \
        cat {input.dict} - | samtools view -F0x900 -b - > {output.bam} 2>> {log}

        samtools sort {output.bam} -T {params.prefix} --threads {threads} -m {params.m} -o {output.sort} 2>>{log}
        samtools index -m {params.m} -@ 1 {output.sort} {output.bai} 2>>{log}
        
        # counts with header
        echo '{wildcards.barcode}' > {output.counts}
        samtools idxstats {output.sort} | grep -v "*" | cut -f3 >> {output.counts}
        """

def get_qout(wildcards, type_o):
    barcodes = get_demux_barcodes(wildcards)
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
    conda: "../envs/minimap2.yaml"
    log: "logs/quant/rowname_kOTU.log"
    benchmark: "benchmarks/quant/rowname_kOTU.txt"
    shell:
        """
        echo '#OTU ID' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule count_matrix:
    input:
        rowname_seqs = rules.rowname_kOTU.output,
        seqs_count = lambda wildcards: get_qout(wildcards, "count"),
    output: "count_matrix.tsv"
    shell:
        """
        cp {input.rowname_seqs} {output}
        # split if exceed 1000 (max 1024 by default in most OS)
        while read fs; do
            paste {output} $fs > {output}.tmp
            # redirection first
            mv {output}.tmp {output}
        done < <(echo {input.seqs_count} | xargs -n 1000)
        """

rule q2_repseqs_import:
    input: rules.rename_fasta_header.output
    output: temp("rep_seqs.qza")
    conda: "../envs/q2plugs.yaml"
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
    conda: "../envs/q2plugs.yaml"
    log: "logs/chimeraF/q2_ftable_import.log"
    benchmark: "benchmarks/chimeraF/q2_ftable_import.txt"
    shell:
        """
        biom convert -i {input} -o {output.biom} --to-hdf5 > {log} 2>&1
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
    conda: "../envs/q2plugs.yaml"
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
    conda: "../envs/q2plugs.yaml"
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
    conda: "../envs/q2plugs.yaml"
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
    conda: "../envs/q2plugs.yaml"
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
    conda: "../envs/q2plugs.yaml"
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