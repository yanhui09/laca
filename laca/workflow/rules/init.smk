check_list_ele("classifier", config["classifier"], ["q2blast", "kraken2", "mmseqs2"])
# rules to initialize workflows, e.g., building database, etc.
DATABASE_DIR = config["database_dir"].rstrip("/")

def get_classifier(c = config["classifier"][0]):
    if c == "kraken2":
        out = expand(DATABASE_DIR + "/kraken2/{prefix}.k2d", prefix = ["hash", "opts", "taxo"])
    elif c == "mmseqs2":
        out = expand(DATABASE_DIR + "/mmseqs2/customDB/customDB{ext}", ext = ["_mapping", "_taxonomy"])
    elif c == "q2blast":
        out = expand(DATABASE_DIR + "/rescript/silva_{prefix}.qza", prefix = ["seqs", "tax"])
    return out 

rule initDB:
    input: get_classifier()

# use predefined MMseqs2 database
rule databases_mmseqs2:
    output: 
        taxdb = expand(DATABASE_DIR + "/mmseqs2/" + config["mmseqs2"]["taxdb"] + "/" + config["mmseqs2"]["taxdb"] + "{ext}",
         ext = ["", ".dbtype", ".index", ".lookup", ".source", ".version",
                "_h", "_h.dbtype", "_h.index", "_mapping", "_taxonomy"]),
        tmp = temp(directory(DATABASE_DIR + "/mmseqs2/tmp")),
    message: "Downloading MMseqs2 database [{params.taxdb}]"
    conda: "../envs/mmseqs2.yaml"
    params:
        taxdb = config["mmseqs2"]["taxdb"],
        targetDB = DATABASE_DIR + "/mmseqs2/" + config["mmseqs2"]["taxdb"] + "/" + config["mmseqs2"]["taxdb"],
    log: "logs/taxonomy/mmseqs2/databases_mmseqs2.log"
    benchmark: "benchmarks/taxonomy/mmseqs2/databases_mmseqs2.txt"
    shell: 
        "mmseqs databases {params.taxdb} {params.targetDB} {output.tmp} 1> {log} 2>&1"

# create seqconfig["mmseqs2"]["taxdb"]manually
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
        blastdb_alias = config["mmseqs2"]["blastdb_alias"],
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
        db = DATABASE_DIR + "/mmseqs2/customDB/blastdb/" + config["mmseqs2"]["blastdb_alias"],
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

# q2 classify-consensus-blast
# rescript silva 138.1, --p-mode majority, --p-perc-identity 0.99
rule rescript_get_silva:
    output: 
        seqs = temp(DATABASE_DIR + "/rescript/silva_rna_seqs.qza"),
        tax = temp(DATABASE_DIR + "/rescript/silva_tax.qza"),
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/rescript/get_silva.log"
    benchmark: "benchmarks/taxonomy/q2blast/rescript/get_silva.txt"
    shell:
        """
        qiime rescript get-silva-data \
            --p-version '138.1' \
            --p-target 'SSURef_NR99' \
            --p-include-species-labels \
            --p-no-rank-propagation \
            --o-silva-sequences {output.seqs} \
            --o-silva-taxonomy {output.tax} 1> {log} 2>&1
        """

rule rescript_reverse_transcribe:
    input: rules.rescript_get_silva.output.seqs
    output: temp(DATABASE_DIR + "/rescript/silva_seqs.qza")
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/rescript/reverse_transcribe.log"
    benchmark: "benchmarks/taxonomy/q2blast/rescript/reverse_transcribe.txt"
    shell:
        """
        qiime rescript reverse-transcribe \
            --i-rna-sequences {input} \
            --o-dna-sequences {output} 1> {log} 2>&1
        """

rule rescript_cull:
    input: rules.rescript_reverse_transcribe.output
    output: temp(DATABASE_DIR + "/rescript/silva_seqs_culled.qza")
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/rescript/cull.log"
    benchmark: "benchmarks/taxonomy/q2blast/rescript/cull.txt"
    threads: config["threads"]["large"]
    shell:
        """
        qiime rescript cull-seqs \
            --i-sequences {input} \
            --o-clean-sequences {output} \
            --p-n-jobs {threads} 1> {log} 2>&1
        """

rule rescript_filter:
    input: 
        seqs = rules.rescript_cull.output,
        tax = rules.rescript_get_silva.output.tax,
    output: 
        filt = temp(DATABASE_DIR + "/rescript/silva_seqs_filt.qza"),
        discard = temp(DATABASE_DIR + "/rescript/silva_seqs_discard.qza")
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/rescript/filter.log"
    benchmark: "benchmarks/taxonomy/q2blast/rescript/filter.txt"
    shell:
        """
        qiime rescript filter-seqs-length-by-taxon \
            --i-sequences {input.seqs} \
            --i-taxonomy {input.tax} \
            --p-labels Archaea Bacteria Eukaryota \
            --p-min-lens 900 1200 1400 \
            --o-filtered-seqs {output.filt} \
            --o-discarded-seqs {output.discard} 1> {log} 2>&1
        """

rule rescript_drep:
    input:
        seqs = rules.rescript_filter.output.filt,
        tax = rules.rescript_get_silva.output.tax,
    output:
        seqs = DATABASE_DIR + "/rescript/silva_seqs_drep.qza",
        tax = DATABASE_DIR + "/rescript/silva_tax_drep.qza",
    conda: "../envs/rescript.yaml"
    log: "logs/taxonomy/q2blast/rescript/drep.log"
    benchmark: "benchmarks/taxonomy/q2blast/rescript/drep.txt"
    threads: config["threads"]["large"]
    shell:
        """
        qiime rescript dereplicate \
            --i-sequences {input.seqs} \
            --i-taxa {input.tax} \
            --p-perc-identity 1.0 \
            --p-rank-handles 'silva' \
            --p-mode 'uniq' \
            --o-dereplicated-sequences {output.seqs} \
            --o-dereplicated-taxa {output.tax} \
            --p-threads {threads} 1> {log} 2>&1
        """