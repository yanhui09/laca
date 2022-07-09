rule guppy:
    input: INPUT_DIR
    output: temp(touch(".demultiplexed"))
    log: "logs/demultiplex_guppy.log"
    benchmark: "benchmarks/demultiplex_guppy.txt"
    threads: config["threads"]["large"]
    params:
        barcode_kits=config["barcode_kits"],
        dir=os.path.join(os.getcwd(), "demultiplexed_guppy"),
    shell: "{workflow.basedir}/resources/ont-guppy/bin/guppy_barcoder -i {input} -s {params.dir} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}"

checkpoint demultiplex_check:
    input: ancient(".demultiplexed")
    output: directory("demultiplexed")
    log: "logs/demultiplex_check.log"
    benchmark: "benchmarks/demultiplex_check.txt"
    params:
        dir=os.path.join(os.getcwd(), "demultiplexed_guppy"),
        nreads_m=config["nreads_m"],
    shell: 
        """
        mv {params.dir} {output} > {log} 2>&1
        # rm shallow sequencing to avoid bardcode bleeding
        mkdir {output}/suspected -p >> {log} 2>&1
        for i in {output}/*/
        do
            if [ "$i" = "{output}/suspected/" ] || [ "$i" = "{output}/unclassified/" ]
            then
                continue
            fi
            nlines=$(cat $i/* | wc -l)
            nreads=$((nlines / 4))
            if [ $nreads -lt {params.nreads_m} ]
            then
                mv "$i" {output}/suspected/ >> {log} 2>&1
                echo "$i moved to {output}/suspected/ due to shallow sequencing < {params.nreads_m}" >> {log} 2>&1
            fi
        done
        """

# collect demultiplexed files
rule collect_fastq:
    input:  "demultiplexed/{barcode}"
    output: temp("qc/{barcode}.fastq")
    shell: "cat {input}/*.fastq > {output}"

def get_demultiplexed(wildcards):
    barcodes = glob_wildcards(checkpoints.demultiplex_check.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    return sorted(set(barcodes))

# consider add minibar? https://github.com/calacademy-research/minibar
# demultiplex? https://github.com/jfjlaros/demultiplex