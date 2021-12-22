# build tree with sepp
# download database
rule download_silva_sepp:
    input: OUTPUT_DIR + "/rep_seq.fasta"
    output: DATABASE_DIR + "/sepp/sepp_silva128.qza"
    message: "Downloading SILVA database for SEPP"
    params:
        link = "https://data.qiime2.org/2021.11/common/sepp-refs-silva-128.qza",
    log: OUTPUT_DIR + "/logs/download_silva_sepp.log"
    benchmark: OUTPUT_DIR + "/benchmarks/download_silva_sepp.txt"
    shell: 
        "wget -O {output} {params.link} 1> {log} 2>&1"

rule q2_rep_seqs:
    input: OUTPUT_DIR + "/rep_seq.fasta"
    output: OUTPUT_DIR + "/sepp/rep_seqs.qza"
    conda: "../envs/q2_sepp.yaml"
    log: OUTPUT_DIR + "/logs/q2_rep_seqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/q2_rep_seqs.txt"
    shell:
        """
        qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path {input} \
        --output-path {output} \
        1> {log} 2>&1
        """

rule q2_sepp:
    input: 
        ref = rules.download_silva_sepp.output,
        rep = rules.q2_rep_seqs.output,
    output:
        tre = OUTPUT_DIR + "/sepp/sepp_tree.qza",
        place = OUTPUT_DIR + "/sepp/sepp_place.qza",
    conda: "../envs/q2_sepp.yaml"
    log: OUTPUT_DIR + "/logs/q2_sepp.log"
    benchmark: OUTPUT_DIR + "/benchmarks/q2_sepp.txt"
    threads: config["threads"]["large"]
    shell:
        """
        qiime fragment-insertion sepp \
        --i-representative-sequences {input.rep} \
        --i-reference-database {input.ref} \
        --o-tree {output.tre} \
        --o-placements {output.place} \
        --p-threads {threads} \
        """

rule q2_ftable:
    input: OUTPUT_DIR + "/count_matrix.tsv"
    output: 
        biom = OUTPUT_DIR + "/sepp/ftable.biom",
        qza = OUTPUT_DIR + "/sepp/ftable.qza",
    conda: "../envs/q2_sepp.yaml"
    log: OUTPUT_DIR + "/logs/q2_ftable.log"
    benchmark: OUTPUT_DIR + "/benchmarks/q2_ftable.txt"
    shell:
        """
        biom convert -i {input} -o {output.biom} --to-hdf5 --table-type="OTU table" 1> {log} 2>&1
        qiime tools import \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --input-path {output.biom} \
        --output-path {output.qza} \
        1>> {log} 2>&1
        """

rule q2_filter:
    input:
        tre = rules.q2_sepp.output.tre,
        table = rules.q2_ftable.output.qza,
    output:
        filtered = OUTPUT_DIR + "/sepp/ftable_filtered.qza",
        removed = OUTPUT_DIR + "/sepp/ftable_removed.qza",
    conda: "../envs/q2_sepp.yaml"
    log: OUTPUT_DIR + "/logs/q2_filter.log"
    benchmark: OUTPUT_DIR + "/benchmarks/q2_filter.txt"
    shell:
       """
       qiime fragment-insertion filter-features \
       --i-tree {input.tre} \
       --i-table {input.table} \
       --o-filtered-table {output.filtered} \
       --o-removed-table {output.removed} \
       """

rule q2_export:
    input:
        filtered = rules.q2_filter.output.filtered,
        tre = rules.q2_sepp.output.tre,
    output:
        filtered_biom = OUTPUT_DIR + "/sepp/ftable_filtered.biom",
        filtered_tsv = OUTPUT_DIR + "/sepp/ftable_filtered.tsv",
        tre = OUTPUT_DIR + "/sepp/sepp.tre",
    conda: "../envs/q2_sepp.yaml"
    log: OUTPUT_DIR + "/logs/q2_export.log"
    benchmark: OUTPUT_DIR + "/benchmarks/q2_export.txt"
    shell:
        """
        qiime tools export \
        --input-path {input.filtered} \
        --output-path {output.filtered_biom} \
        --output-format BIOMV210Format \
        1> {log} 2>&1
        biom convert -i {output.filtered_biom} -o {output.filtered_tsv} --to-tsv 1>> {log} 2>&1
        qiime tools export \
        --input-path {input.tre} \
        --output-path {output.tre} \
        --output-format Newick \
        1>> {log} 2>&1
        """

# only take the aligned sequences
# max fprimer and min rprimer
fprimers_max = config["fprimer_max"]
rprimers_min = config["rprimer_min"]
# pattern
f5al_pattern1 = linked_pattern(fprimers_max, rprimers_min)
f5al_pattern2 = linked_pattern(rprimers_min, fprimers_max)
f5al_patterns = f5al_pattern1 + ' ' + f5al_pattern2

rule trim_repseqs:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: OUTPUT_DIR + "/tree/rep_seqs_trimmed.fasta"
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5al_patterns,
    log: OUTPUT_DIR + "/logs/tree/trim_repseqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/tree/trim_repseqs.txt"
    threads: config["threads"]["normal"]
    shell: "cutadapt -j {threads} {params.f} -o {output} {input} > {log} 2>&1"

rule q2_repseqs:
    input: rules.trim_repseqs.output
    output: OUTPUT_DIR + "/tree/rep_seqs.qza"
    conda: "../envs/q2_phylogen.yaml"
    log: OUTPUT_DIR + "/logs/tree/q2_repseqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/tree/q2_repseqs.txt"
    shell:
        """
        qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path {input} \
        --output-path {output} \
        > {log} 2>&1
        """

rule q2_fasttree:
    input: rules.q2_repseqs.output
    output: 
        expand(OUTPUT_DIR + "/tree/FastTree/{prefix}.qza",
         prefix = ["masked_alignment", "rooted_tree", "alignment", "tree"]),
        _dir = directory(OUTPUT_DIR + "/tree/FastTree"),
    conda: "../envs/q2_phylogen.yaml"
    params:
        supp = config["q2-phylogen"]["fasttree"]
    log: OUTPUT_DIR + "/logs/tree/q2_fasttree.log"
    benchmark: OUTPUT_DIR + "/benchmarks/tree/q2_fasttree.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences {input} \
        --output-dir {output._dir} \
        --p-n-threads {threads} \
        {params.supp} > {log} 2>&1 
        """

rule q2_iqtree:
    input: rules.q2_repseqs.output
    output: 
        expand(OUTPUT_DIR + "/tree/IQ-TREE/{prefix}.qza",
         prefix = ["masked_alignment", "rooted_tree", "alignment", "tree"]),
        _dir = directory(OUTPUT_DIR + "/tree/IQ-TREE"),
    conda: "../envs/q2_phylogen.yaml"
    params:
        supp = config["q2-phylogen"]["iqtree"]
    log: OUTPUT_DIR + "/logs/tree/q2_iqtree.log"
    benchmark: OUTPUT_DIR + "/benchmarks/tree/q2_iqtree.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        qiime phylogeny align-to-tree-iqtree \
        --i-sequences {input} \
        --output-dir {output._dir} \
        --p-n-threads {threads} \
        {params.supp} > {log} 2>&1 
        """

rule q2_raxml:
    input: rules.q2_repseqs.output
    output: 
        expand(OUTPUT_DIR + "/tree/RAxML/{prefix}.qza",
         prefix = ["masked_alignment", "rooted_tree", "alignment", "tree"]),
        _dir = directory(OUTPUT_DIR + "/tree/RAxML"),
    conda: "../envs/q2_phylogen.yaml"
    params:
        supp = config["q2-phylogen"]["raxml"]
    log: OUTPUT_DIR + "/logs/tree/q2_raxml.log"
    benchmark: OUTPUT_DIR + "/benchmarks/tree/q2_raxml.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        qiime phylogeny align-to-tree-raxml \
        --i-sequences {input} \
        --output-dir {output._dir} \
        --p-n-threads {threads} \
        {params.supp} > {log} 2>&1 
        """

rule q2export_tree:
    input: OUTPUT_DIR + "/tree/{phylogen}/rooted_tree.qza"
    output: OUTPUT_DIR + "/tree/{phylogen}/rooted_tree.nwk"
    conda: "../envs/q2_phylogen.yaml"
    log: OUTPUT_DIR + "/logs/tree/{phylogen}_export.log"
    benchmark: OUTPUT_DIR + "/benchmarks/tree/{phylogen}_export.txt"
    shell:
        """
        qiime tools export \
        --input-path {input} \
        --output-path {output} \
        --output-format Newick \
        > {log} 2>&1
        """

rule get_tree:
    input: OUTPUT_DIR + "/tree/" + config["phylogen"][0] + "/rooted_tree.nwk",
    output: OUTPUT_DIR + "/tree.nwk"
    shell:
        "cp -f {input} {output}"