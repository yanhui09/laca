# check phylogen choice
check_list_ele("phylogeny", config["phylogeny"], ["FastTree", "IQ-TREE", "RAxML"])

# only take the aligned sequences
# max fprimer and min rprimer
fprimers_max = config["fprimer_max"]
rprimers_min = config["rprimer_min"]
# pattern
f5al_pattern1 = linked_pattern(fprimers_max, rprimers_min)
f5al_pattern2 = linked_pattern(rprimers_min, fprimers_max)
f5al_patterns = f5al_pattern1 + ' ' + f5al_pattern2

rule trim_repseqs:
    input: get_repseqs()
    output: temp("tree/rep_seqs_trimmed.fasta")
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5al_patterns,
    log: "logs/tree/trim_repseqs.log"
    benchmark: "benchmarks/tree/trim_repseqs.txt"
    threads: config["threads"]["normal"]
    shell: "cutadapt -j {threads} {params.f} -o {output} {input} > {log} 2>&1"

def trim_check2(
    trim = config["trim"], 
    bascdir = config["basecalled_dir"], 
    demuxdir = config["demultiplexed_dir"], 
    merge_runs = config["merge_runs"],
    uchime = config["uchime"]
    ):
    check_val("trim", trim, bool)
    out = rules.trim_repseqs.output
    if trim ==  False:
        out = get_repseqs(bascdir, demuxdir, merge_runs, uchime)
    return out

rule q2_repseqs:
    input: trim_check2()
    output: temp("tree/rep_seqs.qza")
    conda: "../envs/q2plugs.yaml"
    log: "logs/tree/q2_repseqs.log"
    benchmark: "benchmarks/tree/q2_repseqs.txt"
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
        temp(
            expand(
                "tree/FastTree/{prefix}.qza", 
                prefix = ["alignment", "masked_alignment", "tree", "rooted_tree"]
                )
        ),
    conda: "../envs/q2plugs.yaml"
    params:
        supp = config["q2phylo"]["fasttree"]
    log: "logs/tree/q2_fasttree.log"
    benchmark: "benchmarks/tree/q2_fasttree.txt"
    threads: config["threads"]["large"]
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences {input} \
        --o-alignment {output[0]} \
        --o-masked-alignment {output[1]} \
        --o-tree {output[2]} \
        --o-rooted-tree {output[3]} \
        --p-n-threads {threads} \
        {params.supp} > {log} 2>&1 
        """

rule q2_iqtree:
    input: rules.q2_repseqs.output
    output:
        temp(
            expand(
                "tree/IQ-TREE/{prefix}.qza", 
                prefix = ["alignment", "masked_alignment", "tree", "rooted_tree"]
                )
            ),
    conda: "../envs/q2plugs.yaml"
    params:
        supp = config["q2phylo"]["iqtree"]
    log: "logs/tree/q2_iqtree.log"
    benchmark: "benchmarks/tree/q2_iqtree.txt"
    threads: config["threads"]["large"]
    shell:
        """
        qiime phylogeny align-to-tree-mafft-iqtree \
        --i-sequences {input} \
        --o-alignment {output[0]} \
        --o-masked-alignment {output[1]} \
        --o-tree {output[2]} \
        --o-rooted-tree {output[3]} \
        --p-n-threads {threads} \
        {params.supp} > {log} 2>&1 
        """

rule q2_raxml:
    input: rules.q2_repseqs.output
    output:
        temp(
            expand(
                "tree/RAxML/{prefix}.qza", 
                prefix = ["alignment", "masked_alignment", "tree", "rooted_tree"]
                )
            ),
    conda: "../envs/q2plugs.yaml"
    params:
        supp = config["q2phylo"]["raxml"]
    log: "logs/tree/q2_raxml.log"
    benchmark: "benchmarks/tree/q2_raxml.txt"
    threads: config["threads"]["large"]
    shell:
        """
        qiime phylogeny align-to-tree-mafft-raxml \
        --i-sequences {input} \
        --o-alignment {output[0]} \
        --o-masked-alignment {output[1]} \
        --o-tree {output[2]} \
        --o-rooted-tree {output[3]} \
        --p-n-threads {threads} \
        {params.supp} > {log} 2>&1 
        """

rule q2export_tree:
    input: "tree/{phylogen}/rooted_tree.qza"
    output: "tree/{phylogen}/tree.nwk"
    conda: "../envs/q2plugs.yaml"
    params:
        _dir = "tree/{phylogen}",
    log: "logs/tree/{phylogen}_export.log"
    benchmark: "benchmarks/tree/{phylogen}_export.txt"
    shell:
        """
        qiime tools export \
        --input-path {input} \
        --output-path {params._dir} \
        > {log} 2>&1
        """
localrules: get_tree
rule get_tree:
    input: "tree/" + config["phylogeny"][0] + "/tree.nwk",
    output: "tree.nwk"
    shell:
        "cp -f {input} {output}"