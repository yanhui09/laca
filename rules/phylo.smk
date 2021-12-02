# LCA taxonomy with MMseqs2
# download MMseqs2 database SILVA
#mmseqs databases SILVA silvadb tmp
rule download_silva_mmseqs:
    input: rules.count_matrix.output
    output: DATABASE_DIR + "/mmseq2/silvadb"
    message: "Downloading SILVA database for MMseqs2"
    params:
        tmp = DATABASE_DIR + "/mmseq2/tmp",
    conda: "../envs/mmseqs2.yml"
    log: OUTPUT_DIR + "/logs/download_silva.log"
    benchmark: OUTPUT_DIR + "/benchmarks/download_silva.txt"
    shell: 
        "mmseqs databases SILVA {output} {params.tmp} 1> {log} 2>&1"

rule create_queryDB:
    input: rules.dereplicate_denoised_seqs.output
    output: OUTPUT_DIR + "/mmseq2/queryDB"
    conda: "../envs/mmseqs2.yml"
    log: OUTPUT_DIR + "/logs/create_queryDB.log"
    benchmark: OUTPUT_DIR + "/benchmarks/create_queryDB.txt"
    shell: 
        "mmseqs createdb {input} {output} 1> {log} 2>&1"

rule taxonomy_mmseqs2:
    input:
        queryDB = rules.create_queryDB.output,
        targetDB = rules.download_silva_mmseqs.output,
    output:
        taxaDB = OUTPUT_DIR + "/mmseq2/resultDB",
        tsv = OUTPUT_DIR + "/mmseq2/lca_taxonomy.tsv",
    message: "Running LCA taxonomy with MMseqs2"
    params:
        tmp = OUTPUT_DIR + "/mmseq2/tmp",
    conda: "../envs/mmseqs2.yml"
    log: OUTPUT_DIR + "/logs/taxonomy_mmseqs2.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy_mmseqs2.txt"
    threads: config["threads"]["large"]
    shell:
        """
        mmseqs taxonomy {input.queryDB} {input.targetDB} {output.resultDB} {params.tmp} --threads {threads} 1> {log} 2>&1
        mmseqs createtsv {output.resultDB} {output.tsv} 1>> {log} 2>&1
        """

# build tree with seep
# download database
rule download_silva_seep:
    input: rules.count_matrix.output
    output: DATABASE_DIR + "/seep/sepp_silva128.qza"
    message: "Downloading SILVA database for SEEP"
    params:
        link = "https://data.qiime2.org/2021.11/common/sepp-refs-silva-128.qza",
    log: OUTPUT_DIR + "/logs/download_silva.log"
    benchmark: OUTPUT_DIR + "/benchmarks/download_silva.txt"
    shell: 
        "wget -O {output} {params.link} 1> {log} 2>&1"



