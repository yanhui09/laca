# check classifier choice
check_list_ele("classifier", config["classifier"], ["q2blast", "kraken2", "mmseqs2"])

# LCA taxonomy with MMseqs2
rule createdb_query:
    input: get_repseqs()
    output: 
        expand("taxonomy/mmseqs2/queryDB{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source",
                "_h", "_h.dbtype", "_h.index"]),
    conda: "../envs/mmseqs2.yaml"
    params:
        DB = "taxonomy/mmseqs2/queryDB",
    log: "logs/taxonomy/mmseqs2/create_queryDB.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/create_queryDB.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: "mmseqs createdb {input} {params.DB} 1> {log} 2>&1"

def get_seqTaxDB(blastdb_alias, db):
    if blastdb_alias == "":
        dbname = db
    else:
        dbname = "customDB"
    return expand(DATABASE_DIR + "/mmseqs2/"+ dbname + "/" + dbname + "{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source",
                "_h", "_h.dbtype", "_h.index", "_mapping", "_taxonomy"])

rule classify_mmseqs2:
    input:
        ancient(get_seqTaxDB(config["mmseqs2"]["blastdb_alias"], config["mmseqs2"]["taxdb"])[-0]),
        queryDB = rules.createdb_query.output[0],
        targetDB = ancient(get_seqTaxDB(config["mmseqs2"]["blastdb_alias"], config["mmseqs2"]["taxdb"])[0]),
    output:
        resultDB = multiext("taxonomy/mmseqs2/resultDB",
         ".0", ".1", ".2", ".3", ".4", ".5", ".dbtype", ".index"),
        tmp = temp(directory("taxonomy/mmseqs2/tmp")),
    message: "Taxonomy assignment with MMseqs2"
    conda: "../envs/mmseqs2.yaml"
    params: 
        resultDB = "taxonomy/mmseqs2/resultDB",
        lca_mode = config["mmseqs2"]["lca"],
    log: "logs/taxonomy/mmseqs2/taxonomy.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/taxonomy.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        "mmseqs taxonomy {input.queryDB} {input.targetDB} {params.resultDB} {output.tmp}"
        " --search-type 3 --lca-mode {params.lca_mode} --tax-lineage 1"
        " --threads {threads} 1> {log} 2>&1"

rule createtsv:
    input:
        queryDB = rules.createdb_query.output[0],
        resultDB = rules.classify_mmseqs2.output.resultDB,
    output: temp("taxonomy/mmseqs2/output.tsv"),
    message: "Creating tsv file for taxonomy"
    conda: "../envs/mmseqs2.yaml"
    params: 
        resultDB = "taxonomy/mmseqs2/resultDB",
    log: "logs/taxonomy/mmseqs2/createtsv_mmseqs2.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/createtsv_mmseqs2.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        "mmseqs createtsv  {input.queryDB} {params.resultDB} {output}"
        " --full-header --threads {threads}"
        " 1> {log} 2>&1"

localrules: taxonomy_mmseqs2_taxonkit, taxonomy_mmseqs2_silva, taxonomy_mmseqs2, taxonomy_kraken2, repseqs_split, col_q2blast_batch, get_taxonomy
# use taxonkit to get NCBI taxonomy, not for silva
rule lineage_taxonkit_mmseqs2:
    input:
        taxdump = ancient(rules.update_taxdump.output),
        tsv = rules.createtsv.output,
    output: "taxonomy/mmseqs2/lineage.tsv",
    conda: "../envs/taxonkit.yaml"
    params:
        f = 2,
    log: "logs/taxonomy/mmseqs2/lineage_taxonkit.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/lineage_taxonkit.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        "cut -f {params.f} {input.tsv} | sort -u | taxonkit reformat -I 1 -P -j {threads}"
        " --data-dir {input.taxdump} -o {output} 1> {log} 2>&1"

rule taxonomy_mmseqs2_taxonkit:
    input:
        rules.createtsv.output,
        rules.lineage_taxonkit_mmseqs2.output,
    output: "taxonomy/mmseqs2/taxonomy_taxonkit.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep = "\t", header = None)
        df = df.iloc[:, 0:2]
        df.columns = ['kOTUid', 'taxid']
        tax = pd.read_csv(input[1], sep = "\t", header = None)
        tax.columns = ['taxid', 'lineage']
        df = df.merge(tax, how = "left",  on = 'taxid')
        df = df[['kOTUid', 'lineage']]
        df.rename(columns={'kOTUid': 'kOTU_ID', 'lineage': 'taxonomy'}, inplace=True)
        df.to_csv(output[0], sep = "\t", index = False, header = False)

rule taxonomy_mmseqs2_silva:
    input: rules.createtsv.output
    output: "taxonomy/mmseqs2/taxonomy_silva.tsv"
    shell: "cut -f 1,5 {input} | sed 's/\"//g' > {output}"

def get_taxonomy_mmseqs2(taxdb=config["mmseqs2"]["taxdb"], blastdb=config["mmseqs2"]["blastdb_alias"]):
    if taxdb == "silva" and blastdb == "":
        return rules.taxonomy_mmseqs2_silva.output
    else:
        return rules.taxonomy_mmseqs2_taxonkit.output

rule taxonomy_mmseqs2:
    input: get_taxonomy_mmseqs2()
    output: "taxonomy/mmseqs2/taxonomy.tsv"
    shell: "cp -f {input} {output}" 

rule classify_kraken2:
    input:
        ancient(expand(DATABASE_DIR + "/kraken2/{prefix}.k2d", prefix = ["hash", "opts", "taxo"])), 
        fna = ancient(get_repseqs()),
    output: temp("taxonomy/kraken2/classified.tsv"),
    conda: "../envs/kraken2.yaml"
    params:
        classify_cmd = config["kraken2"]["classify_cmd"],
        dbloc = ancient(DATABASE_DIR + "/kraken2"),
    log: "logs/taxonomy/kraken2/classify.log"
    benchmark: "benchmarks/taxonomy/kraken2/classify.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["simple"],
    shell:
        "kraken2 --db {params.dbloc} --threads {threads} --output {output} {input.fna}"
        " {params.classify_cmd} 1> {log} 2>&1"

# use taxonkit to get NCBI taxonomy 
# Not work with special database: greengenes, silva, and RDP
use rule lineage_taxonkit_mmseqs2 as lineage_taxonkit_kraken2 with:
    input:
        taxdump = ancient(rules.update_taxdump.output),
        tsv = rules.classify_kraken2.output,
    output: 
        temp("taxonomy/kraken2/lineage.tsv"),
    params:
        f = 3,
    log: 
        "logs/taxonomy/kraken2/lineage_taxonkit.log"
    benchmark: 
        "benchmarks/taxonomy/kraken2/lineage_taxonkit.txt"

rule taxonomy_kraken2:
    input:
        rules.classify_kraken2.output,
        rules.lineage_taxonkit_kraken2.output,
    output: "taxonomy/kraken2/taxonomy.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep = "\t", header = None)
        df = df.iloc[:, 1:3]
        df.columns = ['kOTUid', 'taxid']
        tax = pd.read_csv(input[1], sep = "\t", header = None)
        tax.columns = ['taxid', 'lineage']
        df = df.merge(tax, how = "left",  on = 'taxid')
        df = df[['kOTUid', 'lineage']]
        df.to_csv(output[0], sep = "\t", index = False, header = False)

# custom parallel, split repseqs into N parts with at most 50 sequnces
checkpoint repseqs_split:
    input: get_repseqs()
    output: temp(directory("taxonomy/q2blast/split"))
    params:
        s = config["q2blast"]["split_by_nseq"],
    log: "logs/taxonomy/q2blast/repseqs_split.log"
    benchmark: "benchmarks/taxonomy/q2blast/repseqs_split.txt"
    shell: "seqkit split {input} -s {params.s} -O {output} 1> {log} 2>&1"

rule q2import:
    input: "taxonomy/q2blast/split/{part}.fasta"
    output: temp("taxonomy/q2blast/split/{part}.qza")
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/q2import_{part}.log"
    benchmark: "benchmarks/taxonomy/q2blast/q2import_{part}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        qiime tools import \
            --input-path {input} \
            --output-path {output} \
            --type 'FeatureData[Sequence]' \
            1> {log} 2>&1
        """
# For NFS disk in clusters
# export TMPDIR=$LOCALSCRATCH in qsub script
#https://forum.qiime2.org/t/error-when-executing-qiime-tools-import-script-on-a-server/7790/2
#https://forum.qiime2.org/t/qiime-tools-import-on-cluster/17531/30

rule classify_q2blast:
    input: 
        query_fna = "taxonomy/q2blast/split/{part}.qza",
        ref_seqs = ancient(rules.rescript_drep.output.seqs),
        ref_tax = ancient(rules.rescript_drep.output.tax),
    output:
        tax = temp("taxonomy/q2blast/blast/{part}_tax.qza"),
        hits = temp("taxonomy/q2blast/blast/{part}_hits.qza"),
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/q2classify_{part}.log"
    benchmark: "benchmarks/taxonomy/q2blast/q2classify_{part}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        qiime feature-classifier classify-consensus-blast \
            --i-query {input.query_fna} \
            --i-reference-reads {input.ref_seqs} \
            --i-reference-taxonomy {input.ref_tax} \
            --o-classification {output.tax} \
            --o-search-results {output.hits} \
            1> {log} 2>&1
        """

rule q2export_tax:
    input: rules.classify_q2blast.output.tax
    output: temp("taxonomy/q2blast/blast/{part}_tax.tsv")
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/q2export_tax_{part}.log"
    benchmark: "benchmarks/taxonomy/q2blast/q2export_tax_{part}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {output} \
            --output-format 'TSVTaxonomyFormat' \
            1> {log} 2>&1
        """    

rule q2export_hits:
    input: rules.classify_q2blast.output.hits
    output: temp("taxonomy/q2blast/blast/{part}_hits.tsv")
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/q2export_hits_{part}.log"
    benchmark: "benchmarks/taxonomy/q2blast/q2export_hits_{part}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {output} \
            --output-format 'BLAST6Format' \
            1> {log} 2>&1
        """

def get_q2blast_batch(wildcards):
    parts = glob_wildcards(checkpoints.repseqs_split.get(**wildcards).output[0] + "/{part}.fasta").part
    return parts

rule col_q2blast_batch:
    input:
        "taxonomy/q2blast/split", 
        tax = lambda wc : expand("taxonomy/q2blast/blast/{part}_tax.tsv", part = get_q2blast_batch(wc)),
        hits = lambda wc : expand("taxonomy/q2blast/blast/{part}_hits.tsv", part = get_q2blast_batch(wc)),
    output:
        tax = "taxonomy/q2blast/taxonomy.tsv",
        hits = "taxonomy/q2blast/hits.tsv",
    shell:
        """
        cat {input.tax} | grep -v "Feature ID" | sort -k1 -V | sed -e 's/; /;/g' -e 's/d__/k__/' > {output.tax}
        cat {input.hits} | sort -k1 -V > {output.hits}
        """

rule get_taxonomy:
    input: "taxonomy/" + config["classifier"][0] + "/taxonomy.tsv",
    output: "taxonomy.tsv"
    shell: "cut -f1,2 {input} > {output}"