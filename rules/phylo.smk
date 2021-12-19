# LCA taxonomy with MMseqs2
# download MMseqs2 database SILVA
DB_NAME = config["mmseqs"]["dbname"]
rule database_mmseqs2:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: 
        db = expand(DATABASE_DIR + "/mmseqs2/"+ DB_NAME + "/" + DB_NAME + "{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source", ".version",
                "_h", "_h.dbtype", "_h.index", "_mapping", "_taxonomy"]),
        tmp = temp(directory(DATABASE_DIR + "/mmseqs2/tmp")),
    message: "Downloading MMseqs2 database [{params.dbname}]"
    conda: "../envs/mmseqs2.yaml"
    params:
        dbname = DB_NAME,
        targetDB = DATABASE_DIR + "/mmseqs2/"+ DB_NAME + "/" + DB_NAME,
    log: OUTPUT_DIR + "/logs/download_silva_mmseqs.log"
    benchmark: OUTPUT_DIR + "/benchmarks/download_silva_mmseqs.txt"
    shell: 
        "mmseqs databases {params.dbname} {params.targetDB} {output.tmp} 1> {log} 2>&1"

rule createdb_mmseqs2:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: 
        expand(OUTPUT_DIR + "/mmseqs2/queryDB{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source",
                "_h", "_h.dbtype", "_h.index"]),
    conda: "../envs/mmseqs2.yaml"
    params:
        queryDB = OUTPUT_DIR + "/mmseqs2/queryDB",
    log: OUTPUT_DIR + "/logs/create_queryDB.log"
    benchmark: OUTPUT_DIR + "/benchmarks/create_queryDB.txt"
    shell: 
        "mmseqs createdb {input} {params.queryDB} 1> {log} 2>&1"

rule taxonomy_mmseqs2:
    input:
        queryDB = rules.createdb_mmseqs2.output[0],
        targetDB = rules.database_mmseqs2.output.db[0],
    output:
        resultDB = multiext(OUTPUT_DIR + "/mmseqs2/resultDB",
         ".0", ".1", ".2", ".3", ".4", ".5", ".dbtype", ".index"),
        tmp = temp(directory(OUTPUT_DIR + "/mmseqs2/tmp")),
    message: "Running LCA taxonomy with MMseqs2"
    conda: "../envs/mmseqs2.yaml"
    params: 
        resultDB = OUTPUT_DIR + "/mmseqs2/resultDB",
        lca_mode = config["mmseqs"]["lca-mode"],
    log: OUTPUT_DIR + "/logs/taxonomy_mmseqs2.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy_mmseqs2.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs taxonomy {input.queryDB} {input.targetDB} {params.resultDB} {output.tmp}"
       # " --lca-ranks species,genus,family,order,class,phylum,superkingdom"
        " --search-type 3 --lca-mode {params.lca_mode} --tax-lineage 1"
        " --threads {threads} 1> {log} 2>&1"

rule createtsv_mmseqs2:
    input:
        queryDB = rules.createdb_mmseqs2.output[0],
        resultDB = rules.taxonomy_mmseqs2.output.resultDB,
    output: OUTPUT_DIR + "/mmseqs2/taxonomy.tsv",
    message: "Creating tsv file for LCA taxonomy"
    conda: "../envs/mmseqs2.yaml"
    params: 
        resultDB = OUTPUT_DIR + "/mmseqs2/resultDB",
    log: OUTPUT_DIR + "/logs/createtsv_mmseqs2.log"
    benchmark: OUTPUT_DIR + "/benchmarks/createtsv_mmseqs2.txt"
    threads: config["threads"]["normal"]
    shell:
        "mmseqs createtsv  {input.queryDB} {params.resultDB} {output}"
        " --full-header --threads {threads}"
        " 1> {log} 2>&1"

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