import os
import re
#--------------
configfile: "config.yaml"
wildcard_constraints:
        barcode="*BRK[0-9][0-9]$"

#--------------
INPUT_DIR = config["basecalled_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
LIVE_BATCH = glob_wildcards(INPUT_DIR + "/{batchid}.fastq").batchid
#---------------

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all:
    input: 
        #expand(OUTPUT_DIR + "/bioms_out/{live_batch}.biom", live_batch = LIVE_BATCH),
        OUTPUT_DIR + "/ONT-L7-GG.txt"

checkpoint demultiplex:
    input:  
        INPUT_DIR + "/{live_batch}.fastq"
    output: 
        temp(directory(OUTPUT_DIR + "/demultiplexed_fq/{live_batch}"))
    threads: 4
    shell:
        """
        mkdir -p {output} && cp {input} {output}
        ./{config[guppy_dir]}/guppy_barcoder -i {output} -s {output} -t {threads} --barcode_kits SQK16S-GXO
        """

rule nanofilt:
    input: 
        OUTPUT_DIR + "/demultiplexed_fq/{live_batch}/{barcode}"
    output: 
        temp(OUTPUT_DIR + "/filt_fq/{live_batch}/{barcode}.fastq")
    conda:
        "workflow/envs/nanofilt.yaml"
    shell: 
        """
        cat {input}/*.fastq | NanoFilt -q 8 -l 1000 --maxlength 1600 --headcrop 15 --tailcrop 15 > {output}
        """

rule fq2fa:
    input: 
        OUTPUT_DIR + "/filt_fq/{live_batch}/{barcode}.fastq"
    output: 
        temp(OUTPUT_DIR + "/filt_fa/{live_batch}/{barcode}.fasta")
    shell: 
        """
        config[usearch11] -fastq_filter {input} -fastaout {output}
        # find filt_fa/{wildcards.live_batch}/*.fasta -size -100k -delete # remove small files 0-10 reads
        """ 

rule fa_adjust:
    input: 
        OUTPUT_DIR + "/filt_fa/{live_batch}/{barcode}.fasta"
    output: 
        temp(OUTPUT_DIR + "/filt_fa/{live_batch}/adjusted-{barcode}.fasta")
    shell: 
        "scripts/fasta_number.py {input} {wildcards.barcode}_ > {output}"

rule qiime1_assignment:
    input: 
        OUTPUT_DIR + "/filt_fa/{live_batch}/adjusted-{barcode}.fasta"
    output:  
        temp(directory(OUTPUT_DIR + "/uclust/{live_batch}/{barcode}"))
    conda: 
        "envs/qiime1.yaml"
    shell: 
        """
        if [ ! -f config[gg_13_8_dir]/rep_set/99_otus.fasta ]; then; 
            gzip -d config[gg_13_8_dir]/rep_set/99_otus.fasta.gz
        fi
        
        if [ ! -f config[gg_13_8_dir]/taxonomy/99_otu_taxonomy.txt ]; then; 
            gzip -d config[gg_13_8_dir]/taxonomy/99_otu_taxonomy.txt.gz
        fi
        
        parallel_assign_taxonomy_uclust.py -i {input} -r config[gg_13_8_dir]/rep_set/99_otus.fasta \
        -t config[gg_13_8_dir]/taxonomy/99_otu_taxonomy.txt \
        -o {output} -O 1 --similarity 0.80 --min_consensus_fraction 0.60
        """

rule biom_per_barcode:
    input: 
        OUTPUT_DIR + "/uclust/{live_batch}/{barcode}"
    output: 
        temp(OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}.biom") 
    conda:
        "envs/qiime1.yaml"
    shell:
        """
        echo {wildcards.live_batch}
        echo {wildcards.barcode}
        
        awk 'BEGIN {FS = "\t"}; {print $2}'  {input}/adjusted-{wildcards.barcode}_tax_assignments.txt > {input}/taxa.tmp
        num=$(wc -l taxa.tmp | awk '{print $1}')-1
        for ((i=0; i<=$num;i++)); do echo 1; done > {input}/col2
        for ((i=0; i<=$num;i++)); do echo ONTU_$i; done > {input}/col1
        
        paste {input}/col1 {input}/col2 {input}/taxa.tmp > {input}/step1.tmp
        
        TAB=$'\t'
        HASH=$'#'
        echo "$HASH OTU_ID $TAB {wildcards.barcode} $TAB taxonomy" | sed "s/ //g" > {input}/tmp.txt
        cat {input}/table.tmp {input}/step1.tmp | sed "s/^$/Unassigned/g"> {input}/{wildcards.barcode}.tmp
        
        biom convert -i {input}/{wildcards.barcode}.tmp -o {output}  --process-obs-metadata=taxonomy --table-type="OTU table" --to-json 
        """

rule taxa_summary:
    input: 
        OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}.biom"
    output: 
        temp(OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}_L7.biom")
    params: 
        out_dir = OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}_L7"
    conda:
        "envs/qiime1.yaml"    
    shell: 
        """
        summarize_taxa.py -i {input} -L 7 -a -o {params.out_dir} -u 0.0005  # remove bellow 0.05 % as most of it is usually mess
        mv {params.out_dir}/*.biom {output}
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output[0]
    return expand(OUTPUT_DIR + "/live_bioms/{live_batch}/{barcode}_L7.biom",
           live_batch = wildcards.live_batch, 
           barcode = glob_wildcards(checkpoint_output + "/{barcode}/{runid}.fastq").barcode)

rule biom_per_batch:
    input: 
        aggregate_input
    output: 
        OUTPUT_DIR + "/bioms_out/{live_batch}.biom"
    conda: 
        "envs/qiime1.yaml"
    shell:
        "merge_otu_tables.py -i {input} -o {output}"

rule biom_update:
    input: 
        expand(OUTPUT_DIR + "/bioms_out/{live_batch}.biom", live_batch = LIVE_BATCH)
    output: 
        biom_l7 = OUTPUT_DIR + "/ONT-L7-GG.biom",
        tsv_l7 = OUTPUT_DIR + "/ONT-L7-GG.txt"
    conda:
        "envs/qiime1.yaml" 
    shell:
        """
        merge_otu_tables.py -i {input} -o {output.biom_l7}
        biom convert -i {output.biom_l7} -o {output.tsv_l7} --to-tsv
        """