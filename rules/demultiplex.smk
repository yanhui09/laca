checkpoint guppy_demultiplex:
    input: INPUT_DIR
    output: directory(OUTPUT_DIR + "/demultiplexed")
    log: OUTPUT_DIR + "/logs/demultiplex.log"
    threads: config["threads"]["large"]
    params:
        guppy=config["guppy"],
        barcode_kits=config["barcode_kits"],
    shell:
        "{params.guppy}/guppy_barcoder -i {input} -s {output} -t {threads} --barcode_kits {params.barcode_kits} 2>{log}"

# collect demultiplexed files
rule collect_fastq:
    input:  OUTPUT_DIR + "/demultiplexed/{barcode}"
    output: temp(OUTPUT_DIR + "/raw/{barcode}.fastq")
    shell: "cat {input}/*.fastq > {output}"

def get_demultiplexed(wildcards):
    barcodes = glob_wildcards(checkpoints.guppy_demultiplex.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    return list(set(barcodes))

# consider add minibar? https://github.com/calacademy-research/minibar