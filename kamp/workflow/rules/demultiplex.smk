rule guppy:
    input: INPUT_DIR
    output: ".demultiplexed"
    log: "logs/demultiplex.log"
    threads: config["threads"]["large"]
    params:
        guppy=config["guppy"],
        barcode_kits=config["barcode_kits"],
        out="demultiplexed",
    shell:
        "{params.guppy}/guppy_barcoder -i {input} -s {params.out} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}"

checkpoint demultiplex_check:
    input: ".demultiplexed"
    output: touch(directory("demultiplexed"))

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