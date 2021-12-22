# check phylogen choice
def check_phylogen_val(phylogen):
    if phylogen:
        for i in phylogen:
            if i not in ["FastTree", "IQ-TREE", "RAxML"]:
                raise ValueError("\t\nPhylogenetic method to generate the tree file not recognized.\t\nPlease choose from FastTree, IQ-TREE, or RXaML in the config.yaml file.")
    else:
        raise ValueError("\t\nPhylogenetic method to generate the tree file not specified.\t\nPlease choose from FastTree, IQ-TREE, or RXaML in the config.yaml file.")

check_phylogen_val(config["phylogen"])

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