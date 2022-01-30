# pseudo checkpoint for isONcorCon
checkpoint clusters_isONcor:
    input: OUTPUT_DIR + "/isONclustCon/clusters"
    output: directory(OUTPUT_DIR + "/isONcorCon/clusters")
    shell: "cp -rf {input} {output}"

# run_isoncorrect provide multithread processing for isONcorrect in batches
# racon seems not to be supported in batch mode
rule isONcorrect:
    input: rules.get_fqs_split.output,
    output:
        _dir = directory(OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/{clust_id}/isONcor"), 
        fq = OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/{clust_id}/isONcor/{clust_id}/corrected_reads.fastq"
    conda: "../envs/isONcorCon.yaml"
    params:
        _dir = OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/{clust_id}",
    log: OUTPUT_DIR + "/logs/isONcorCon/{barcode}/{c}/{clust_id}/isONcorrect.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorCon/{barcode}/{c}/{clust_id}/isONcorrect.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        mkdir -p {params._dir}/isONclust
        cp {input} {params._dir}/isONclust
        run_isoncorrect --t {threads} --fastq_folder {params._dir}/isONclust --outfolder {output._dir} \
        --set_w_dynamically --split_wrt_batches > {log} 2>&1
        rm -rf {params._dir}/isONclust
        """

rule isoCon:
    input: rules.isONcorrect.output.fq
    output:
        _dir = directory(OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/{clust_id}/isoCon"),
        fna = OUTPUT_DIR + "/isONcorCon/{barcode}/{c}/{clust_id}/isoCon/final_candidates.fa",
    conda: "../envs/isONcorCon.yaml"
    log: OUTPUT_DIR + "/logs/isONcorCon/{barcode}/{c}/{clust_id}/isoCon.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONcorCon/{barcode}/{c}/{clust_id}/isoCon.txt"
    threads: config["threads"]["normal"]
    shell: "IsoCon pipeline -fl_reads {input} -outfolder {output._dir} --nr_cores {threads} > {log} 2>&1"
    
def get_isONcorCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.clusters_isONcor.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    fnas = []
    for i in bc_kb_cis:
        bc, kb, ci = i.split("_")
        fnas.append(OUTPUT_DIR + "/isONcorCon/{bc}/{kb}/{ci}/isoCon/final_candidates.fa".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_isONcorCon:
    input: lambda wc: get_isONcorCon(wc)
    output: OUTPUT_DIR + "/isONcorCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i, id_i = [ i.split("/")[index] for index in [-5, -4, -3] ]
                with open(i, "r") as inp:
                    j = 0
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "_cand" + str(j) + "\n"
                            j += 1
                        out.write(line)
