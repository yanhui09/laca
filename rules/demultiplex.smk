rule guppy:
    input: INPUT_DIR
    output: OUTPUT_DIR + "/.demultiplexed"
    log: OUTPUT_DIR + "/logs/demultiplex.log"
    threads: config["threads"]["large"]
    params:
        guppy=config["guppy"],
        barcode_kits=config["barcode_kits"],
        out=OUTPUT_DIR + "/demultiplexed",
    shell:
        "{params.guppy}/guppy_barcoder -i {input} -s {params.out} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}"

checkpoint demultiplex_check:
    input: OUTPUT_DIR + "/.demultiplexed"
    output: touch(directory(OUTPUT_DIR + "/demultiplexed"))

# collect demultiplexed files
rule collect_fastq:
    input:  OUTPUT_DIR + "/demultiplexed/{barcode}"
    output: temp(OUTPUT_DIR + "/qc/{barcode}.fastq")
    shell: "cat {input}/*.fastq > {output}"

def get_demultiplexed(wildcards):
    barcodes = glob_wildcards(checkpoints.demultiplex_check.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    return sorted(set(barcodes))

# consider add minibar? https://github.com/calacademy-research/minibar