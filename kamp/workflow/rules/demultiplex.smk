rule guppy:
    input: INPUT_DIR
    output: directory("demultiplexed_guppy")
    log: "logs/demultiplex_guppy.log"
    benchmark: "benchmarks/demultiplex_guppy.txt"
    threads: config["threads"]["large"]
    params:
        barcode_kits=config["barcode_kits"],
    shell: "{workflow.basedir}/resources/ont-guppy/bin/guppy_barcoder -i {input} -s {output} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}"

checkpoint demultiplex_check:
    input: directory("demultiplexed_guppy")
    output: directory("demultiplexed")
    shell: "mv {input} {output}"

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