rule isONclust:
    input: rules.split_by_cluster.output
    output:
        _dir = directory(OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/isONclust"),
        tsv = OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/isONclust/final_clusters.tsv"
    conda: "../envs/isONcorCon.yaml"
    log: OUTPUT_DIR + "/logs/isONcorCon/{barcode}/{c}/isONclust/cluster.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorCon/{barcode}/{c}/isONclust/cluster.txt"
    threads: config["threads"]["normal"]
    shell:
        "isONclust --ont --fastq {input} "
        "--outfolder {output._dir} --t {threads} > {log} 2>&1"

rule isONclust_write:
    input:
        tsv = rules.isONclust.output.tsv,
        fqs = rules.split_by_cluster.output,
    output: directory(OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/isONclust/fastq_files")
    conda: "../envs/isONcorCon.yaml"
    log: OUTPUT_DIR + "/logs/isONcorCon/{barcode}/{c}/isONclust/write.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorCon/{barcode}/{c}/isONclust/write.txt"
    shell:
        "isONclust write_fastq --N 1 "
        "--clusters {input.tsv} --fastq {input.fqs} "
        "--outfolder  {output} > {log} 2>&1"

# rm snakemake_timestamp as the run_isoncorrect take ids through str.split(".")
# https://github.com/ksahlin/isONcorrect/blob/master/run_isoncorrect#L153
rule clean_timestamp:
    input: rules.isONclust_write.output
    output: OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/isONclust/.timestamp_removed"
    shell:
        """
        rm -f {input}/.snakemake_timestamp
        touch {output}
        """

# run_isoncorrect provide multithread processing for isONcorrect in batches
# racon seems not to be supported in batch mode
checkpoint isONcorrect:
    input:
        rules.clean_timestamp.output,
        fqs = rules.isONclust_write.output,
    output: directory(OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/isONcorrect"),
    conda: "../envs/isONcorCon.yaml"
    log: OUTPUT_DIR + "/logs/isONcorCon/{barcode}/{c}/isONcorrect.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorCon/{barcode}/{c}/isONcorrect.txt"
    threads: config["threads"]["normal"]
    shell:
        "run_isoncorrect --t {threads} --fastq_folder {input.fqs} --outfolder {output} "
        "--set_w_dynamically --split_wrt_batches "
        "> {log} 2>&1"

rule IsoCon:
    input: OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/isONcorrect/{cid}/corrected_reads.fastq"
    output:
        _dir = directory(OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/IsoCon/id_{cid}"),
        fna = OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/IsoCon/id_{cid}/final_candidates.fa",
    conda: "../envs/isONcorCon.yaml"
    log: OUTPUT_DIR + "/logs/isONcorCon/{barcode}/{c}/IsoCon/id_{cid}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorCon/{barcode}/{c}/IsoCon/id_{cid}.txt"
    threads: config["threads"]["normal"]
    shell:
        "IsoCon pipeline -fl_reads {input} -outfolder {output._dir} --nr_cores {threads} > {log} 2>&1"
    
def get_isONcorCon(wildcards, pooling = True):
    check_val("pooling", pooling, bool)
    if pooling == True:
        barcodes = ["pooled"]
    else:
        barcodes = get_demultiplexed(wildcards)

    fnas = []
    for i in barcodes:
        cs = glob_wildcards(checkpoints.cluster_info.get(barcode=i).output[0] + "/{c}.txt").c
        for j in cs:
            cids = glob_wildcards(checkpoints.isONcorrect.get(barcode=i, c=j).output[0] + "/{cid}/corrected_reads.fastq").cid
            for k in cids:
                fnas.append(OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/IsoCon/id_{cid}/final_candidates.fa".format(barcode=i, c=j, cid=k))
    return fnas

rule collect_isONcorCon:
    input: lambda wc: get_isONcorCon(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/isONcorCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i, id_i = [ i.split("/")[index] for index in [-5, -4, -2] ]
                with open(i, "r") as inp:
                    j = 0
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "_cand" + str(j) + "\n"
                            j += 1
                        out.write(line)
