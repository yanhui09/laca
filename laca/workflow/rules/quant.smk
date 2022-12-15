# check classifier choice
check_list_ele("cluster", config["cluster"], ["kmerCon", "clustCon", "NGSpeciesID", "NGSpeciesID2", "isONcorCon", "umiCon"])
# check quantification method
check_list_ele("quant", config["quant"], ["seqid", "minimap2"])

def get_consensus2derep(cls):
   #{cls}/{cls}.fna except umiCon ({cls}/{cls}_trimmed.fna)
    if cls == "umiCon":
        return "{cls}/{cls}_trimmed.fna".format(cls=cls)
    else:
        return "{cls}/{cls}.fna".format(cls=cls)

# dereplicate sequences with MMseqs2
rule drep_consensus:
    input: [get_consensus2derep(cls=i) for i in config["cluster"]]
    output: 
        rep = temp("quant/mmseqs2_rep_seq.fasta"),
        all_by_cluster = temp("quant/mmseqs2_all_seqs.fasta"),
        tsv = temp("quant/mmseqs2_cluster.tsv"),
        tmp = temp(directory("quant/tmp")),
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
        "mmseqs easy-cluster {input[0]} {params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.mid} -c {params.c} --cluster-reassign > {log} 2>&1"

# keep fasta header unique
rule rename_drep_seqs:
    input: rules.drep_consensus.output.rep
    output: "rep_seqs.fasta"
    run:
        with open(output[0], "w") as out:
            with open (input[0], "r") as inp:
                i = 1
                for line in inp:
                    if line.startswith(">"):
                        line = ">OTU_" + str(i) + "\n"
                        i += 1 
                    out.write(line)

def get_cand_cls(cls=config["cluster"][0]):
    if cls == "kmerCon":
        return "kmerBin/clusters"
    else:
        return "{cls}/clusters".format(cls=cls)

rule combine_cls:
    input: get_cand_cls()
    output: temp("quant/cls_combined.tsv")
    run:
        import pandas as pd
        # append cand suffix by cls
        cls = input[0].split("/")[-2]
        if cls == "kmerBin":
            cand = "_0cand1"
        elif cls == "isONclustCon" or cls == "clustCon" or cls == "umiCon":
            cand = "cand1"
        else:
            cand = ""

        # combine csv files under input
        files = os.listdir(input[0])
        cls_csvs = [os.path.join(input[0], i) for i in files if i.endswith(".csv")]
        # concatanate all csv files, with file name (without suffix) as column
        for i in cls_csvs:
            df = pd.read_csv(i, sep="\t", header=None)
            df["cls"] = os.path.basename(i).split(".")[0] + cand
            if i == cls_csvs[0]:
                df_all = df
            else:
                df_all = pd.concat([df_all, df], axis=0)
        # write to output
        df_all.to_csv(output[0], sep="\t", header=None, index=None)

# abudance matrix by read id
rule matrix_seqid: 
    input:
        rules.drep_consensus.output.tsv,
        rules.combine_cls.output,
        fqs = lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_qced_barcodes(wc)),
    output: "quant/matrix_seqid.tsv"
    run:
        import pandas as pd
        # OTU <- derep_cls -> cand_cls
        derep_cls = pd.read_csv(input[0], sep="\t", header=None)
        derep_cls.columns = ["rep_cls","cls"]
        # OTU with incremental number by rep_cls, 1, 2, ...
        derep_cls["OTU"] = derep_cls["rep_cls"].factorize()[0] + 1

        cand_cls = pd.read_csv(input[1], sep="\t", header=None)
        cand_cls.columns = ["seqid", "cls"]
        # merge on cls, left join, dummy files not included in {cls}/{cls}.fna 
        cls = derep_cls.merge(cand_cls, on="cls", how="left")
        # take header from fastq files, ^@, as a list
        for fq in input.fqs:
            barcode = os.path.basename(fq).split(".")[0]
            seqid = []
            with open(fq, "r") as f:
                for line in f:
                    if line.startswith("@"):
                        seqid.append(line.rstrip().split(" ")[0][1:])
            # if cls["seqid"] is in seqid, then 1, else 0
            cls[barcode] = cls["seqid"].isin(seqid).astype(int)
        
        # exclude rep_cls, cls, seqid
        # group by OTU, sum by barcode
        cls = cls.drop(["rep_cls", "cls", "seqid"], axis=1).groupby("OTU").sum()
        # append 'OTU_' after sorting
        cls.index = ["OTU_" + str(i) for i in cls.index.sort_values()]
        # rename OTU as "# OTU ID" for qiime2 import
        cls.index.name = "#OTU ID"
        # write to output
        cls.to_csv(output[0], sep="\t", header=True, index=True)
                
# abundance matrix with minimap2
rule index_repseqs:
    input: rules.rename_drep_seqs.output
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
        max_de = config["minimap2"]["max_de"],
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
rule rowname_OTU:
    input:
        bam = lambda wildcards: get_qout(wildcards, "bam"),
        bai = lambda wildcards: get_qout(wildcards, "bai"),
    output: temp("quant/rowname_seqs")
    conda: "../envs/minimap2.yaml"
    log: "logs/quant/rowname_OTU.log"
    benchmark: "benchmarks/quant/rowname_OTU.txt"
    shell:
        """
        echo '#OTU ID' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule matrix_minimap2:
    input:
        rowname_seqs = rules.rowname_OTU.output,
        seqs_count = lambda wildcards: get_qout(wildcards, "count"),
    output: "quant/matrix_minimap2.tsv"
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

rule count_matrix:
    input: ["quant/matrix_{quant}.tsv".format(quant=i) for i in config["quant"]], 
    output: "count_matrix.tsv"
    shell: "cp {input[0]} {output}"

rule q2repseqs_import:
    input: rules.rename_drep_seqs.output
    output: temp("rep_seqs.qza")
    conda: "../envs/q2plugs.yaml"
    log: "logs/uchime/q2_repseqs_import.log"
    benchmark: "benchmarks/uchime/q2_repseqs_import.txt"
    shell:
        """
        qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path {input} \
        --output-path {output} \
        > {log} 2>&1
        """

rule q2ftable_import:
    input: rules.count_matrix.output
    output: 
        biom = temp("table.biom"),
        qza = temp("table.qza")
    conda: "../envs/q2plugs.yaml"
    log: "logs/uchime/q2_ftable_import.log"
    benchmark: "benchmarks/uchime/q2_ftable_import.txt"
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

rule q2uchime_denovo:
    input:
        ftable = rules.q2ftable_import.output.qza,
        repseqs = rules.q2repseqs_import.output,
    output: 
        nonchimeras = "uchime/nonchimeras.qza",
        chimeras = "uchime/chimeras.qza",
        stats = "uchime/stats.qza",
    conda: "../envs/q2plugs.yaml"
    log: "logs/uchime/uchime_denovo.log"
    benchmark: "benchmarks/uchime/uchime_denovo.txt"
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

rule q2filter_features:
    input:
        ftable = rules.q2ftable_import.output.qza,
        nonchimeras = rules.q2uchime_denovo.output.nonchimeras,
    output: temp("uchime/table_nonchimeras.qza")
    conda: "../envs/q2plugs.yaml"
    log: "logs/uchime/filter_features.log"
    benchmark: "benchmarks/uchime/filter_features.txt"
    shell:
        """
        qiime feature-table filter-features \
        --i-table {input.ftable} \
        --m-metadata-file {input.nonchimeras} \
        --o-filtered-table {output} \
        > {log} 2>&1
        """

rule q2filter_seqs:
    input:
        repseqs = rules.q2repseqs_import.output,
        nonchimeras = rules.q2uchime_denovo.output.nonchimeras,
    output: temp("uchime/rep_seqs_nonchimeras.qza")
    conda: "../envs/q2plugs.yaml"
    log: "logs/uchime/filter_seqs.log"
    benchmark: "benchmarks/uchime/filter_seqs.txt"
    shell:
        """
        qiime feature-table filter-seqs \
        --i-data {input.repseqs} \
        --m-metadata-file {input.nonchimeras} \
        --o-filtered-data {output} \
        > {log} 2>&1
        """

rule q2ftable_export:
    input: rules.q2filter_features.output
    output:
        biom = temp("uchime/table_nonchimeras.biom"), 
        tsv = "count_matrix_filt.tsv"
    params:
        _dir = "uchime"
    conda: "../envs/q2plugs.yaml"
    log: "logs/uchime/ftable_export.log"
    benchmark: "benchmarks/uchime/ftable_export.txt"
    shell:
        """
        qiime tools export \
        --input-path {input} \
        --output-path {params._dir} \
        > {log} 2>&1
        mv {params._dir}/feature-table.biom {output.biom}
        biom convert -i {output.biom} -o {output.tsv} --to-tsv >> {log} 2>&1
        """

rule q2repseqs_export:
    input: rules.q2filter_seqs.output
    output: "rep_seqs_filt.fasta"
    conda: "../envs/q2plugs.yaml"
    params:
        _dir = "uchime"
    log: "logs/uchime/repseqs_export.log"
    benchmark: "benchmarks/uchime/repseqs_export.txt"
    shell:
        """
        qiime tools export \
        --input-path {input} \
        --output-path {params._dir} \
        > {log} 2>&1
        mv {params._dir}/dna-sequences.fasta {output}
        """

def get_repseqs(
    bascdir = config["basecalled_dir"], 
    demuxdir = config["demultiplexed_dir"], 
    merge_runs = config["merge_runs"],
    uchime = config["uchime"]
    ):
    # rep_seqs_merged.fasta if 'bascdir' and 'demuxdir' are none
    # and 'merge_runs' is not empty
    if bascdir is None and demuxdir is None and merge_runs is not None:
        return "rep_seqs_merged.fasta"
    else:
        return chimera_filt(uchime)[1]