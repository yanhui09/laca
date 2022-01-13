# check classifier choice
check_list_ele("classifier", config["classifier"], ["kraken2", "mmseqs2"])

# LCA taxonomy with MMseqs2
TaxDB = config["mmseqs"]["taxdb"]
# use predefined MMseqs2 database
rule databases_mmseqs2:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: 
        taxdb = expand(DATABASE_DIR + "/mmseqs2/"+ TaxDB + "/" + TaxDB + "{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source", ".version",
                "_h", "_h.dbtype", "_h.index", "_mapping", "_taxonomy"]),
        tmp = temp(directory(DATABASE_DIR + "/mmseqs2/tmp")),
    message: "Downloading MMseqs2 database [{params.taxdb}]"
    conda: "../envs/mmseqs2.yaml"
    params:
        taxdb = TaxDB,
        targetDB = DATABASE_DIR + "/mmseqs2/"+ TaxDB + "/" + TaxDB,
    log: OUTPUT_DIR + "/logs/taxonomy/mmseqs2/databases_mmseqs2.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/mmseqs2/databases_mmseqs2.txt"
    shell: 
        "mmseqs databases {params.taxdb} {params.targetDB} {output.tmp} 1> {log} 2>&1"

# create seqTaxDB manually
# the predefined one do not have taxonomy information
rule update_taxdump:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: directory(DATABASE_DIR + "/mmseqs2/customDB/taxdump"),
    log: OUTPUT_DIR + "/logs/taxonomy/update_taxdump.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/update_taxdump.txt"
    shell:
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P {output} 1> {log} 2>&1
        tar -xzvf {output}/taxdump.tar.gz -C {output} 1>> {log} 2>&1
        """

rule update_blastdb:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: directory(DATABASE_DIR + "/mmseqs2/customDB/blastdb"),
    conda: "../envs/rsync.yaml"
    params:
        blastdb_alias = config["mmseqs"]["blastdb_alias"],
    log: OUTPUT_DIR + "/logs/taxonomy/update_blastdb.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/update_blastdb.txt"
    shell:
        """
        rsync -rv --include="{params.blastdb_alias}*.tar.gz" --exclude="*" rsync://ftp.ncbi.nlm.nih.gov/blast/db/ {output} 1> {log} 2>&1
        for file in {output}/{params.blastdb_alias}*.tar.gz; do tar -xzvf $file -C {output} 1>> {log} 2>&1; done
        """

rule blastdbcmd:
    input: rules.update_blastdb.output
    output:
        fna = DATABASE_DIR + "/mmseqs2/customDB/custom.fna",
        taxidmapping = DATABASE_DIR + "/mmseqs2/customDB/custom.taxidmapping",
    conda: "../envs/blast.yaml"
    params:
        db = DATABASE_DIR + "/mmseqs2/customDB/blastdb/" + config["mmseqs"]["blastdb_alias"],
    log: OUTPUT_DIR + "/logs/taxonomy/blast/blastdbcmd.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/blast/blastdbcmd.txt"
    shell:
    # ignore err: [blastdbcmd] error while reading seqid
    # I don't know how to escape it using blastdbcmd
        """
        blastdbcmd -db {params.db} -entry all > {output.fna} 2> {log} || true
        blastdbcmd -db {params.db} -entry all -outfmt "%a %T" > {output.taxidmapping} 2>> {log} || true
        """

rule createdb_seqTax:
    input: rules.blastdbcmd.output.fna
    output: 
        expand(DATABASE_DIR + "/mmseqs2/customDB/customDB{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source",
                "_h", "_h.dbtype", "_h.index"]),
    conda: "../envs/mmseqs2.yaml"
    params:
        DB = DATABASE_DIR + "/mmseqs2/customDB/customDB",
    log: OUTPUT_DIR + "/logs/taxonomy/mmseqs2/create_customDB.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/mmseqs2/create_customDB.txt"
    shell: 
        "mmseqs createdb {input} {params.DB} 1> {log} 2>&1"

rule createtaxdb_seqTax:
    input:
        rules.createdb_seqTax.output,
        taxdump = rules.update_taxdump.output,
        taxidmapping = rules.blastdbcmd.output.taxidmapping,
    output: 
        taxdb = expand(DATABASE_DIR + "/mmseqs2/customDB/customDB{ext}",
                     ext = ["_mapping", "_taxonomy"]),
        tmp = temp(directory(DATABASE_DIR + "/mmseqs2/tmp")),
    conda: "../envs/mmseqs2.yaml"
    params:
        DB = DATABASE_DIR + "/mmseqs2/customDB/customDB",
    log: OUTPUT_DIR + "/logs/taxonomy/mmseqs2/create_taxdb.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/mmseqs2/create_taxdb.txt"
    shell: 
        "mmseqs createtaxdb {params.DB} {output.tmp} "
        "--ncbi-tax-dump {input.taxdump} --tax-mapping-file {input.taxidmapping} 1> {log} 2>&1"
     
rule createdb_query:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: 
        expand(OUTPUT_DIR + "/taxonomy/mmseqs2/queryDB{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source",
                "_h", "_h.dbtype", "_h.index"]),
    conda: "../envs/mmseqs2.yaml"
    params:
        DB = OUTPUT_DIR + "/taxonomy/mmseqs2/queryDB",
    log: OUTPUT_DIR + "/logs/taxonomy/mmseqs2/create_queryDB.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/mmseqs2/create_queryDB.txt"
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
        resultDB = multiext(OUTPUT_DIR + "/taxonomy/mmseqs2/resultDB",
         ".0", ".1", ".2", ".3", ".4", ".5", ".dbtype", ".index"),
        tmp = temp(directory(OUTPUT_DIR + "/taxonomy/mmseqs2/tmp")),
    message: "Taxonomy assignment with MMseqs2"
    conda: "../envs/mmseqs2.yaml"
    params: 
        resultDB = OUTPUT_DIR + "/taxonomy/mmseqs2/resultDB",
        lca_mode = config["mmseqs"]["lca-mode"],
    log: OUTPUT_DIR + "/logs/taxonomy/mmseqs2/taxonomy.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/mmseqs2/taxonomy.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs taxonomy {input.queryDB} {input.targetDB} {params.resultDB} {output.tmp}"
        " --search-type 3 --lca-mode {params.lca_mode} --tax-lineage 1"
        " --threads {threads} 1> {log} 2>&1"

rule createtsv:
    input:
        queryDB = rules.createdb_query.output[0],
        resultDB = rules.classify_mmseqs2.output.resultDB,
    output: temp(OUTPUT_DIR + "/taxonomy/mmseqs2/output.tsv"),
    message: "Creating tsv file for taxonomy"
    conda: "../envs/mmseqs2.yaml"
    params: 
        resultDB = OUTPUT_DIR + "/taxonomy/mmseqs2/resultDB",
    log: OUTPUT_DIR + "/logs/taxonomy/mmseqs2/createtsv_mmseqs2.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/mmseqs2/createtsv_mmseqs2.txt"
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
    output: OUTPUT_DIR + "/taxonomy/mmseqs2/lineage.tsv",
    conda: "../envs/taxonkit.yaml"
    params:
        f = 2,
    log: OUTPUT_DIR + "/logs/taxonomy/mmseqs2/lineage_taxonkit.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/mmseqs2/lineage_taxonkit.txt"
    threads: config["threads"]["normal"]
    shell:
        "cut -f {params.f} {input.tsv} | sort -u | taxonkit reformat -I 1 -P -j {threads}"
        " --data-dir {input.taxdump} -o {output} 1> {log} 2>&1"

rule taxonomy_mmseqs2:
    input:
        rules.createtsv.output,
        rules.lineage_taxonkit_mmseqs2.output,
    output: OUTPUT_DIR + "/taxonomy/mmseqs2/taxonomy.tsv"
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

# add kraken classifier
rule build_database:
    input: OUTPUT_DIR + "/rep_seqs.fasta"
    output: 
        #temp(directory(DATABASE_DIR + "/kraken2/library")),
        #temp(directory(DATABASE_DIR + "/kraken2/taxonomy")),
        #temp(DATABASE_DIR + "/kraken2/seqid2taxid.map"),
        k2d = expand(DATABASE_DIR + "/kraken2/{prefix}.k2d", prefix = ["hash", "opts", "taxo"]),
        dbloc = directory(DATABASE_DIR + "/kraken2"),
    conda: "../envs/kraken2.yaml"
    params:
        buildb_cmd = config["kraken2"]["buildb_cmd"],
    log: OUTPUT_DIR + "/logs/taxonomy/kraken2/build_database.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/kraken2/build_database.txt"
    threads: config["threads"]["large"]
    shell:
        "kraken2-build --db {output.dbloc} {params.buildb_cmd} -threads {threads} 1> {log} 2>&1"

rule classify_kraken2:
    input:
        ancient(rules.build_database.output.k2d), 
        dbloc = ancient(DATABASE_DIR + "/kraken2"),
        fna = OUTPUT_DIR + "/rep_seqs.fasta",
    output: OUTPUT_DIR  + "/taxonomy/kraken2/classified.tsv",
    conda: "../envs/kraken2.yaml"
    params:
        classify_cmd = config["kraken2"]["classify_cmd"],
    log: OUTPUT_DIR + "/logs/taxonomy/kraken2/classify.log"
    benchmark: OUTPUT_DIR + "/benchmarks/taxonomy/kraken2/classify.txt"
    threads: config["threads"]["large"]
    shell:
        "kraken2 --db {input.dbloc} --threads {threads} --output {output} {input.fna}"
        " {params.classify_cmd} 1> {log} 2>&1"

# use taxonkit to get NCBI taxonomy
use rule lineage_taxonkit_mmseqs2 as lineage_taxonkit_kraken2 with:
    input:
        taxdump = ancient(rules.update_taxdump.output),
        tsv = rules.classify_kraken2.output,
    output: 
        OUTPUT_DIR + "/taxonomy/kraken2/lineage.tsv",
    params:
        f = 3,
    log: 
        OUTPUT_DIR + "/logs/taxonomy/kraken2/lineage_taxonkit.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/taxonomy/kraken2/lineage_taxonkit.txt"

rule taxonomy_kraken2:
    input:
        rules.classify_kraken2.output,
        rules.lineage_taxonkit_kraken2.output,
    output: OUTPUT_DIR + "/taxonomy/kraken2/taxonomy.tsv"
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
    input: OUTPUT_DIR + "/taxonomy/" + config["classifier"][0] + "/taxonomy.tsv",
    output: OUTPUT_DIR + "/taxonomy.tsv"
    shell:
        "cp -f {input} {output}"