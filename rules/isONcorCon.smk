rule isONclust:
    input: rules.split_by_cluster.output
    output:
        _dir = directory(OUTPUT_DIR + "/isONclust/{barcode}/{c}"),
        tsv = OUTPUT_DIR + "/isONclust/{barcode}/{c}/final_clusters.tsv"
    conda: "../envs/isONcorrect_IsoCon.yaml"
    log: OUTPUT_DIR + "/logs/isONclust/{barcode}/{c}/cluster.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclust/{barcode}/{c}/cluster.txt"
    threads: config["threads"]["normal"]
    shell:
        "isONclust --ont --fastq {input} "
        "--outfolder {output._dir} --t {threads} > {log} 2>&1"

rule isONclust_write:
    input:
        tsv = rules.isONclust.output.tsv,
        fqs = rules.split_by_cluster.output,
    output: directory(OUTPUT_DIR + "/isONclust/{barcode}/{c}/fastq_files")
    conda: "../envs/isONcorrect_IsoCon.yaml"
    log: OUTPUT_DIR + "/logs/isONclust/{barcode}/{c}/write.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclust/{barcode}/{c}/write.txt"
    shell:
        "isONclust write_fastq --N 1 "
        "--clusters {input.tsv} --fastq {input.fqs} "
        "--outfolder  {output} > {log} 2>&1"

checkpoint isONcorrect:
    input: rules.isONclust_write.output,
    output: directory(OUTPUT_DIR + "/isONcorrect/{barcode}/{c}/corrected"),
    conda: "../envs/isONcorrect_IsoCon.yaml"
    log: OUTPUT_DIR + "/logs/isONcorrect/{barcode}/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorrect/{barcode}/{c}.txt"
    threads: config["threads"]["large"]
    shell:
        "run_isoncorrect --t {threads} --fastq_folder {input}  --outfolder {output} > {log} 2>&1"

rule IsoCon:
    input: OUTPUT_DIR + "/isONcorrect/{barcode}/{c}/corrected/{cid}.fastq"
    output:
        _dir = directory(OUTPUT_DIR + "/isONcorrect/{barcode}/{c}/IsoCon/id_{cid}"),
        fna = OUTPUT_DIR + "/IsoCon/{barcode}/{c}/IsoCon/id{cid}/final_candidates.fa"
    conda: "../envs/isONcorrect_IsoCon.yaml"
    log: OUTPUT_DIR + "/logs/isONcorrect/{barcode}/{c}/IsoCon/id_{cid}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorrect/{barcode}/{c}/IsoCon/id_{cid}.txt"
    threads: config["threads"]["large"]
    shell:
        "IsoCon pipeline -fl_reads {input} -outfolder {output} --nr_cores {threads} > {log} 2>&1"
    
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
            cids = glob_wildcards(checkpoints.isONcorrect.get(barcode=i, c=j).output[0] + "/{cid}.fastq").cid
            for k in cids:
                fnas.append(OUTPUT_DIR + "/isONcorrect/{barcode}/{c}/IsoCon/id_{cid}/final_candidates.fa".format(barcode=i, c=j, cid=k))
    return fnas

rule collect_isONcorCon:
    input: lambda wc: get_isONcorCon(wc, pooling = config["pooling"]),
    output: OUTPUT_DIR + "/isONcorCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i, id_i = [ i.split("/")[index] for index in [-5, -4, -2] ]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "\n"
                        out.write(line)
