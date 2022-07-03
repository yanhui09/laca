rule guppy:
    input: INPUT_DIR
    output: touch(".demultiplexed")
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
    params:
        dir=os.path.join(os.getcwd(), "demultiplexed_guppy"),
    shell: "mv {params.dir} {output}"

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