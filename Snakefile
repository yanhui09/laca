# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
INPUT_DIR = config["basecalled_dir"]
OUTPUT_DIR = config["results_dir"]
wildcard_constraints:
        barcode="^BRK"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

#rule all:
#    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.


#include: "rules/common.smk"
#include: "rules/other.smk"
rule demultiplex:
    input: 
        INPUT_DIR + "{live_batch}.fastq"
    output: 
        temp(OUTPUT_DIR + "demultiplexed_fq/{live_batch}")
    threads: 2
    shell:
        "{config[guppy_dir]}/guppy_barcoder -i {input} -s {output}  -t {threads}  --barcode_kits SQK16S-GXO"

rule nanofilt:
    input: 
        OUTPUT_DIR + "demultiplexed_fq/{live_batch}/{barcode}" 
    output: 
        temp(OUTPUT_DIR + "filt_fq/{live_batch}/{barcode}.fastq")
    conda:
        "workflow/envs/nanofilt.yaml"
    shell: 
        """
        cat {input}/*.fastq > {input}.fastq
        cat {input}.fastq | NanoFilt -q 8 -l 1000 --maxlength 1600 --headcrop 15 --tailcrop 15 > {output}
        """

rule fq2fa:
    input: 
        OUTPUT_DIR + "filt_fq/{live_batch}/{barcode}.fastq"
    output: 
        temp(OUTPUT_DIR + "filt_fa/{live_batch}/{barcode}.fasta")
    shell: 
        """
        config[usearch11] -fastq_filter {input} -fastaout {output}
        find filt_fa/{wildcards.live_batch}/*.fasta -size -100k -delete # remove small files 0-10 reads
        """ 

rule fa_number:
    input: 
        OUTPUT_DIR + "filt_fa/{live_batch}/{barcode}.fasta"
    output: 
        temp(OUTPUT_DIR + "filt_fa/{live_batch}/adjusted-{barcode}.fasta")
    conda:
        "envs/fa_number.yaml"
    shell: 
        "scripts/fasta_number.py {input} {wildcards.barcode}_ > {output}

rule qiime1_assignment:
    input: 
        OUTPUT_DIR + "filt_fa/{live_batch}/adjusted-{barcode}.fasta"
    output:  
        temp(OUTPUT_DIR + "uclust/{live_batch}/{barcode}")
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
rule summary:
    input: 
        OUTPUT_DIR + "uclust/{live_batch}/{barcode}"
    output: 
        OUTPUT_DIR + "live_bioms/{live_batch}.biom" 
    conda:
        "envs/qiime1.yaml"
    shell:
        """
        for file in adjusted-_*.txt
        do
        string=${file}
        string2=$(sed "-es/adjusted-filt_//; s/.fastq_tax_assignments.txt//" <<< $string)
        echo $string2
        
        awk 'BEGIN { FS = "\t" }; {print $2}'  $file > taxa.txt
        num=$(wc -l taxa.txt | awk '{print $1}')-1
        for ((i=0; i<=$num;i++)); do echo 1; done > col2
        for ((i=0; i<=$num;i++)); do echo ONTU_$i; done > col1
        paste col1 col2 taxa.txt > step1.txt
        
        rm col1
        rm col2
        rm taxa.txt
        
        TAB=$'\t'
        HASH=$'#'
        echo "$HASH OTU_ID $TAB ${string2} $TAB taxonomy" | sed "s/ //g" > tmp.txt
        cat tmp.txt step1.txt | sed "s/^$/Unassigned/g"> ${string2}.txt
        rm step1.txt
        rm tmp.txt
        
        biom convert -i ${string2}.txt -o ${string2}.biom --process-obs-metadata=taxonomy --table-type="OTU table" --to-json 
        summarize_taxa.py -i *.biom -L 7 -a -o L7 -u 0.0005  #remove bellow 0.05 % as most of it is usually mess
        
        cd ..
        done
        """ 