checkpoint demultiplex:
    input: INPUT_DIR
    output: directory(OUTPUT_DIR + "/demultiplexed")
    log: OUTPUT_DIR + "/logs/demultiplex.log"
    threads: config["threads"]["large"]
    params:
        guppy=config["guppy"],
        barcode_kits=config["barcode_kits"],
        use_cuda=config["use_cuda"],
    shell:
        """
        {params.guppy}/guppy_barcoder -i {input} -s {output} -t {threads} --barcode_kits {params.barcode_kits} 2>{log}
        if {params.use_cuda}; then
            {params.guppy}/guppy_barcoder -i {input} -s {output} -t {threads} --barcode_kits {params.barcode_kits} -x auto 2>{log}
        fi
        """

# collect demultiplexed files
rule collect_fastq:
    input:  OUTPUT_DIR + "/demultiplexed/{barcode}"
    output: temp(OUTPUT_DIR + "/raw_fq/{barcode}.fastq")
    shell: "cat {input}/*.fastq > {output}"
