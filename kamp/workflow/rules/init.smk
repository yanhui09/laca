check_list_ele("classifier", config["classifier"], ["kraken2", "mmseqs2"])

# rules to initialize workflows, e.g., building database, etc.
DATABASE_DIR = config["database_dir"].rstrip("/")

# LCA taxonomy with MMseqs2
TaxDB = config["mmseqs"]["taxdb"]

def get_classifier(c):
    if c == "kraken2":
        out = expand(DATABASE_DIR + "/kraken2/{prefix}.k2d", prefix = ["hash", "opts", "taxo"])
    elif c == "mmseqs2":
        out = expand(DATABASE_DIR + "/mmseqs2/customDB/customDB{ext}", ext = ["_mapping", "_taxonomy"])
    return out 

rule initDB:
    input: get_classifier(config["classifier"][0])
    output: touch(DATABASE_DIR + "/.initDB_DONE")

# use predefined MMseqs2 database
rule databases_mmseqs2:
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
    log: "logs/taxonomy/mmseqs2/databases_mmseqs2.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/databases_mmseqs2.txt"
    shell: 
        "mmseqs databases {params.taxdb} {params.targetDB} {output.tmp} 1> {log} 2>&1"

# create seqTaxDB manually
# the predefined one do not have taxonomy information
rule update_taxdump:
    output: directory(DATABASE_DIR + "/mmseqs2/customDB/taxdump"),
    log: "logs/taxonomy/update_taxdump.log"
    benchmark: "benchmarks/taxonomy/update_taxdump.txt"
    shell:
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P {output} 1> {log} 2>&1
        tar -xzvf {output}/taxdump.tar.gz -C {output} 1>> {log} 2>&1
        """

rule update_blastdb:
    output: directory(DATABASE_DIR + "/mmseqs2/customDB/blastdb"),
    conda: "../envs/rsync.yaml"
    params:
        blastdb_alias = config["mmseqs"]["blastdb_alias"],
    log: "logs/taxonomy/update_blastdb.log"
    benchmark: "benchmarks/taxonomy/update_blastdb.txt"
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
    log: "logs/taxonomy/blast/blastdbcmd.log"
    benchmark: "benchmarks/taxonomy/blast/blastdbcmd.txt"
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
    log: "logs/taxonomy/mmseqs2/create_customDB.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/create_customDB.txt"
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
    log: "logs/taxonomy/mmseqs2/create_taxdb.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/create_taxdb.txt"
    shell: 
        "mmseqs createtaxdb {params.DB} {output.tmp} "
        "--ncbi-tax-dump {input.taxdump} --tax-mapping-file {input.taxidmapping} 1> {log} 2>&1"

# add kraken classifier
#rule database_kraken2:
#    output: 
#        k2d = expand(DATABASE_DIR + "/kraken2/{prefix}.k2d", prefix = ["hash", "opts", "taxo"]),
#    conda: "../envs/kraken2.yaml"
#    params:
#        buildb_cmd = config["kraken2"]["buildb_cmd"],
#        dbloc = directory(DATABASE_DIR + "/kraken2"),
#    log: "logs/taxonomy/kraken2/build_database.log"
#    benchmark: "benchmarks/taxonomy/kraken2/build_database.txt"
#    threads: config["threads"]["large"]
#    shell:
#        "kraken2-build --db {params.dbloc} {params.buildb_cmd} -threads {threads} 1> {log} 2>&1"

# Download prebuilt kraken2 database
rule kraken2_prebuilt:
    output: 
        k2d = expand(DATABASE_DIR + "/kraken2/{prefix}.k2d", prefix = ["hash", "opts", "taxo"]),
    params:
        address = config["kraken2"]["prebuilt"],
        dbloc = directory(DATABASE_DIR + "/kraken2"),
    log: "logs/taxonomy/kraken2/kraken2_prebuilt.log"
    benchmark: "benchmarks/taxonomy/kraken2/kraken2_prebuilt.txt"
    threads: 1
    shell:
        """
        wget {params.address} -O {params.dbloc}/kraken2.tar.gz 1> {log} 2>&1
        tar -xzvf {params.dbloc}/kraken2.tar.gz -C {params.dbloc} 1>> {log} 2>&1
        rm {params.dbloc}/kraken2.tar.gz 1>> {log} 2>&1
        """

