# check classifier choice
check_list_ele("classifier", config["classifier"], ["kraken2", "mmseqs2"])

# LCA taxonomy with MMseqs2
TaxDB = config["mmseqs"]["taxdb"]
    
rule createdb_query:
    input: chimeraF(config["chimeraF"])[1]
    output: 
        expand("taxonomy/mmseqs2/queryDB{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source",
                "_h", "_h.dbtype", "_h.index"]),
    conda: "../envs/mmseqs2.yaml"
    params:
        DB = "taxonomy/mmseqs2/queryDB",
    log: "logs/taxonomy/mmseqs2/create_queryDB.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/create_queryDB.txt"
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
        ancient(get_seqTaxDB(config["mmseqs"]["blastdb_alias"], config["mmseqs"]["taxdb"])[-0]),
        queryDB = rules.createdb_query.output[0],
        targetDB = ancient(get_seqTaxDB(config["mmseqs"]["blastdb_alias"], config["mmseqs"]["taxdb"])[0]),
    output:
        resultDB = multiext("taxonomy/mmseqs2/resultDB",
         ".0", ".1", ".2", ".3", ".4", ".5", ".dbtype", ".index"),
        tmp = temp(directory("taxonomy/mmseqs2/tmp")),
    message: "Taxonomy assignment with MMseqs2"
    conda: "../envs/mmseqs2.yaml"
    params: 
        resultDB = "taxonomy/mmseqs2/resultDB",
        lca_mode = config["mmseqs"]["lca-mode"],
    log: "logs/taxonomy/mmseqs2/taxonomy.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/taxonomy.txt"
    threads: config["threads"]["large"]
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
    shell:
        "mmseqs createtsv  {input.queryDB} {params.resultDB} {output}"
        " --full-header --threads {threads}"
        " 1> {log} 2>&1"

# use taxonkit to get NCBI taxonomy
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
    shell:
        "cut -f {params.f} {input.tsv} | sort -u | taxonkit reformat -I 1 -P -j {threads}"
        " --data-dir {input.taxdump} -o {output} 1> {log} 2>&1"

rule taxonomy_mmseqs2:
    input:
        rules.createtsv.output,
        rules.lineage_taxonkit_mmseqs2.output,
    output: "taxonomy/mmseqs2/taxonomy.tsv"
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

rule classify_kraken2:
    input:
        ancient(DATABASE_DIR + "/.initDB_DONE"),
        ancient(expand(DATABASE_DIR + "/kraken2/{prefix}.k2d", prefix = ["hash", "opts", "taxo"])), 
        fna = chimeraF(config["chimeraF"])[1]
    output: temp("taxonomy/kraken2/classified.tsv"),
    conda: "../envs/kraken2.yaml"
    params:
        classify_cmd = config["kraken2"]["classify_cmd"],
        dbloc = ancient(DATABASE_DIR + "/kraken2"),
    log: "logs/taxonomy/kraken2/classify.log"
    benchmark: "benchmarks/taxonomy/kraken2/classify.txt"
    threads: config["threads"]["large"]
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

rule get_taxonomy:
    input: "taxonomy/" + config["classifier"][0] + "/taxonomy.tsv",
    output: "taxonomy.tsv"
    shell:
        "cp -f {input} {output}"