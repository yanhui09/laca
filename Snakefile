#--------------
configfile: "config.yaml"
wildcard_constraints:
        barcode="BRK[0-9][0-9]"

#--------------
INPUT_DIR = config["basecalled_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
LIVE_BATCH = glob_wildcards(INPUT_DIR + "/{batchid}.fastq").batchid
#---------------

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all:
    input: 
        OUTPUT_DIR + "/ONT-L7-GG.txt"
    
checkpoint demultiplex:
    input:  
        INPUT_DIR + "/{live_batch}.fastq"
    output: 
        temp(directory(OUTPUT_DIR + "/demultiplexed_fq/{live_batch}"))
    log: 
        OUTPUT_DIR + "/logs/demultiplex/{live_batch}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/demultiplex/{live_batch}.tsv"
    threads: 2
    shell:
        """
        mkdir -p {output} && cp {input} {output}
        {config[guppy_dir]}/guppy_barcoder -i {output} -s {output} -t {threads} --barcode_kits SQK16S-GXO -x auto
        """

rule nanofilt:
    input: 
        OUTPUT_DIR + "/demultiplexed_fq/{live_batch}/{barcode}"
    output: 
        temp(OUTPUT_DIR + "/filt_fq/{live_batch}/{barcode}.fastq")
    log: 
        OUTPUT_DIR + "/logs/nanofilt/{live_batch}/{barcode}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/nanofilt/{live_batch}/{barcode}.tsv"
    conda:
        "envs/nanofilt.yaml"
    shell: 
        """
        cat {input}/*.fastq | NanoFilt -q 8 -l 1000 --maxlength 1600 --headcrop 15 --tailcrop 15 > {output}
        """

rule fq2fa:
    input: 
        OUTPUT_DIR + "/filt_fq/{live_batch}/{barcode}.fastq"
    output: 
        temp(OUTPUT_DIR + "/filt_fa/{live_batch}/{barcode}.fasta")
    log: 
        OUTPUT_DIR + "/logs/fq2fa/{live_batch}/{barcode}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/fq2fa/{live_batch}/{barcode}.tsv"
    shell: 
        """
        {config[usearch11]} -fastq_filter {input} -fastaout {output}
        # find filt_fa/{wildcards.live_batch}/*.fasta -size -100k -delete # remove small files 0-10 reads
        """ 

rule fa_adjust:
    input: 
        OUTPUT_DIR + "/filt_fa/{live_batch}/{barcode}.fasta"
    output: 
        temp(OUTPUT_DIR + "/filt_fa/{live_batch}/adjusted-{barcode}.fasta")
    log: 
        OUTPUT_DIR + "/logs/fa_adjust/{live_batch}/{barcode}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/fa_adjust/{live_batch}/{barcode}.tsv"
    shell: 
        "scripts/fasta_number.py {input} {wildcards.barcode}_ > {output}"

rule taxa_assignment:
    input: 
        OUTPUT_DIR + "/filt_fa/{live_batch}/adjusted-{barcode}.fasta"
    output:  
        temp(directory(OUTPUT_DIR + "/uclust/{live_batch}/{barcode}"))
    log: 
        OUTPUT_DIR + "/logs/taxa_assignment/{live_batch}/{barcode}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/taxa_assignment/{live_batch}/{barcode}.tsv"
    threads: 4
    conda: 
        "envs/qiime1.yaml"
    shell: 
        """
        if [ ! -f {config[gg_13_8_dir]}/rep_set/99_otus.fasta ]; then 
            gzip -d {config[gg_13_8_dir]}/rep_set/99_otus.fasta.gz
        fi
        
        if [ ! -f {config[gg_13_8_dir]}/taxonomy/99_otu_taxonomy.txt ]; then 
            gzip -d {config[gg_13_8_dir]}/taxonomy/99_otu_taxonomy.txt.gz
        fi
        
        parallel_assign_taxonomy_uclust.py -i {input} -r {config[gg_13_8_dir]}/rep_set/99_otus.fasta \
        -t {config[gg_13_8_dir]}/taxonomy/99_otu_taxonomy.txt \
        -o {output} -O {threads} --similarity 0.80 --min_consensus_fraction 0.60
        """

rule biom_per_barcode:
    input: 
        OUTPUT_DIR + "/uclust/{live_batch}/{barcode}"
    output: 
        temp(OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}.biom") 
    log: 
        OUTPUT_DIR + "/logs/biom_perbarcode/{live_batch}/{barcode}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/biom_perbarcode/{live_batch}/{barcode}.tsv"
    conda:
        "envs/qiime1.yaml"
    shell:
        """
        scripts/read_count.sh {input} {wildcards.barcode}
        biom convert -i {input}/{wildcards.barcode}.tmp -o {output} --process-obs-metadata=taxonomy --table-type="OTU table" --to-json 
        rm -f {input}/{wildcards.barcode}.tmp 
        """

rule taxa_summary:
    input: 
        OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}.biom"
    output: 
        temp(OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}_L7.biom")
    params:
        dir_out = OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}"
    log: 
        OUTPUT_DIR + "/logs/taxa_summary/{live_batch}/{barcode}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/taxa_summary/{live_batch}/{barcode}.tsv"
    conda:
        "envs/qiime1.yaml"    
    shell: 
        """
        summarize_taxa.py -i {input} -L 7 -a -o {params.dir_out}
        mv {params.dir_out}/{wildcards.barcode}_L7.biom {output}
        rm -rf {params.dir_out}
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output[0]
    return expand(OUTPUT_DIR + "/live_bioms/{{live_batch}}/{barcode}_L7.biom",
           barcode = glob_wildcards(checkpoint_output + "/{barcode, BRK[0-9][0-9]}/{runid}.fastq").barcode)

rule biom_per_batch:
    input: 
        aggregate_input
    output: 
        OUTPUT_DIR + "/bioms_out/{live_batch}.biom"
    params:
        files = lambda wildcards, input: ','.join(input)
    log: 
        OUTPUT_DIR + "/logs/biom_per_batch/{live_batch}.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/biom_per_batch/{live_batch}.tsv"
    conda: 
        "envs/qiime1.yaml"
    shell:
        "merge_otu_tables.py -i {params.files} -o {output}"

rule biom_update:
    input: 
        expand(OUTPUT_DIR + "/bioms_out/{live_batch}.biom", live_batch = LIVE_BATCH)
    output: 
        biom_l7 = OUTPUT_DIR + "/ONT-L7-GG.biom",
        tsv_l7 = OUTPUT_DIR + "/ONT-L7-GG.txt"
    params:
        files = lambda wildcards, input: ','.join(input),
        demux_dir = OUTPUT_DIR + "/demultiplex" 
    log: 
        OUTPUT_DIR + "/logs/biom_update.log"
    #benchmark:
    #    OUTPUT_DIR + "/benchmarks/biom_update.tsv"
    conda:
        "envs/qiime1.yaml" 
    shell:
        """
        merge_otu_tables.py -i {params.files} -o {output.biom_l7}
        biom convert -i {output.biom_l7} -o {output.tsv_l7} --to-tsv
        rm -rf {params.demux_dir} 
        """