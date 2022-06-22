def get_fq4Con(kmerbin = True):
    check_val("kmerbin", kmerbin, bool)
    if kmerbin == True:
        out = rules.split_bin.output
    else:
        out = rules.skip_bin.output
    return out

# kmerCon
rule get_fqs_split:
    input: get_fq4Con(config["kmerbin"]),
    output: OUTPUT_DIR + "/kmerCon/{barcode}/{c}/split/0.fastq",
    log: OUTPUT_DIR + "/logs/kmerCon/{barcode}/{c}/0/get_fqs_split.log"
    benchmark: OUTPUT_DIR + "/benchmarks/kmerCon/{barcode}/{c}/0/get_fqs_split.txt"
    shell: "cp -f {input} {output} 2> {log}"

rule spoa:
    input: rules.get_fqs_split.output
    output: OUTPUT_DIR + "/kmerCon/{barcode}/{c}/polish/0/draft/raw.fna"
    conda: '../envs/spoa.yaml'
    params:
        l = config["spoa"]["l"],
        r = config["spoa"]["r"],
        g = config["spoa"]["g"],
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/0/spoa_consensus.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/0/spoa_consensus.txt"
    shell: "spoa {input} -l {params.l} -r {params.r} -g {params.g} > {output} 2> {log}"

# split kmerbins
# clustCon
# draw draft with max average score from pairwise alignments
rule minimap2clust:
    input: get_fq4Con(config["kmerbin"])
    output: OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/minimap2clust.paf"
    conda: '../envs/polish.yaml'
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/minimap2clust/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/minimap2clust/{c}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x ava-ont --no-long-join -r100"
        " {input} {input} > {output} 2> {log}"

rule bin2clustering:
    input: OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/minimap2clust.paf"
    output: OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust.csv",
    conda: '../envs/clustCon.yaml'
    params:
        prefix = OUTPUT_DIR + "/clustCon/{barcode}/avr_aln/{c}/bin2clust",
        min_score_frac = config["clustCon"]["min_score_frac"],
        min_reads = config["min_cluster_size"],
        max_recurs = config["clustCon"]["max_recursion"],
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/bin2clust/{c}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/bin2clust/{c}.txt"
    shell:
        "python scripts/binClust.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input} > {log} 2>& 1"

def get_clust(wildcards, pool = True, kmerbin = True):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
        
    if kmerbin is True:
        bin2cls = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            bin2cls.append(OUTPUT_DIR + "/clustCon/{bc}/avr_aln/{kb}/bin2clust.csv".format(bc=bc, kb=kb))
    else:
        if pool is True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        bin2cls = expand(OUTPUT_DIR + "/clustCon/{bc}/avr_aln/all/bin2clust.csv", bc=bcs)
    return bin2cls
 
checkpoint cls_clustCon:
    input: lambda wc: get_clust(wc, pool = config["pool"], kmerbin = config["kmerbin"])
    output: directory(OUTPUT_DIR + "/clustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import shutil
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            num_lines = sum(1 for line in open(i))
            if num_lines > params.min_size:
                barcode, c = [ i.split("/")[index] for index in [-4, -2] ]
                df = pd.read_csv(i)
                for clust_id, df_clust in df.groupby('cluster'):
                    if len(df_clust) >= params.min_size:
                        bc_kb_ci = "/{barcode}_{c}_{clust_id}".format(barcode=barcode, c=c, clust_id=clust_id)
                        df_clust['read_id'].to_csv(output[0] + bc_kb_ci + ".csv", header = False, index = False)
                        ref_idx  = df_clust['clust_read_score'].idxmax()
                        ref_read = df_clust.loc[ref_idx, ['read_id']]
                        ref_read.to_csv(output[0] + bc_kb_ci + ".ref", header=False, index=False)
   
rule get_fqs_split2:
    input:
        pool = OUTPUT_DIR + "/clustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        ref = OUTPUT_DIR + "/clustCon/clusters/{barcode}_{c}_{clust_id}.ref",
        binned = get_fq4Con(config["kmerbin"]),
    output:
        pool = OUTPUT_DIR + "/clustCon/{barcode}/{c}/split/{clust_id}.fastq",
        ref = OUTPUT_DIR + "/clustCon/{barcode}/{c}/polish/{clust_id}/draft/raw.fna",
    conda: '../envs/seqkit.yaml'
    log: OUTPUT_DIR + "/logs/clustCon/{barcode}/{c}/{clust_id}/get_fqs_split.log"
    benchmark: OUTPUT_DIR + "/benchmarks/clustCon/{barcode}/{c}/{clust_id}/get_fqs_split.txt"
    shell:
        """
        seqkit grep -f {input.pool} {input.binned} -o {output.pool} --quiet 2> {log}
        seqkit grep -f {input.ref} {input.binned} --quiet | seqkit fq2fa | \
        seqkit replace -p '^(.+)$' -r '{wildcards.barcode}_{wildcards.c}_{wildcards.clust_id}' -o {output.ref} 2> {log}
        """

# isONclust
rule isONclust:
    input: get_fq4Con(config["kmerbin"])
    output:
        _dir = directory(OUTPUT_DIR + "/isONclustCon/{barcode}/{c}"),
        tsv = OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/final_clusters.tsv"
    conda: "../envs/isONcorCon.yaml"
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/isONclust.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/isONclust.txt"
    threads: config["threads"]["large"]
    shell:
        "isONclust --ont --fastq {input} "
        "--outfolder {output._dir} --t {threads} > {log} 2>&1"

def get_isONclust(wildcards, pool = True, kmerbin = True):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
        
    if kmerbin is True:
        bin2cls = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            bin2cls.append(OUTPUT_DIR + "/isONclustCon/{bc}/{kb}/final_clusters.tsv".format(bc=bc, kb=kb))
    else:
        if pool is True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        bin2cls = expand(OUTPUT_DIR + "/isONclustCon/{bc}/all/final_clusters.tsv", bc=bcs)
    return bin2cls

checkpoint cls_isONclust:
    input: lambda wc: get_isONclust(wc, pool = config["pool"], kmerbin = config["kmerbin"])
    output: directory(OUTPUT_DIR + "/isONclustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            barcode, c = [ i.split('/')[index] for index in [-3, -2] ]
            df_i = pd.read_csv(i, sep = '\t', header = None)
            df_i.columns = ['cluster', 'seqid']
            for clust_id, df_clust in df_i.groupby('cluster'):
                if len(df_clust) >= params.min_size:
                    df_clust['seqid'].to_csv(output[0] + "/{barcode}_{c}_{clust_id}.csv".format(barcode=barcode, c=c, clust_id=clust_id),
                     header = False, index = False)

rule get_fqs_split3:
    input:
        bin2clust = OUTPUT_DIR + "/isONclustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        binned = get_fq4Con(config["kmerbin"]),
    output: OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/split/{clust_id}.fastq",
    conda: '../envs/seqkit.yaml'
    log: OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/{clust_id}/get_fqs_split.log"
    benchmark: OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/get_fqs_split.txt"
    shell: "seqkit grep -f {input.bin2clust} {input.binned} -o {output} --quiet 2> {log}"

use rule spoa as spoa2 with:
    input:
        rules.get_fqs_split3.output
    output: 
        OUTPUT_DIR + "/isONclustCon/{barcode}/{c}/polish/{clust_id}/draft/raw.fna"
    log: 
        OUTPUT_DIR + "/logs/isONclustCon/{barcode}/{c}/{clust_id}/spoa.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/spoa.txt"

# polish with racon and medaka
# reused in racon iterations
rule minimap2polish:
    input: 
      ref = OUTPUT_DIR + "/{cls}/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.fna",
      fastq = OUTPUT_DIR + "/{cls}/{barcode}/{c}/split/{clust_id}.fastq",
    output: OUTPUT_DIR + "/{cls}/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.paf",
    message: "Polish draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}]: alignments against {wildcards.assembly} assembly [{wildcards.cls}]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/{cls}/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/{cls}/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.txt"
    threads: config["threads"]["normal"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = OUTPUT_DIR + "/{cls}/{barcode}/{c}/polish/{clust_id}/draft/raw"
        return(prefix + ".paf", prefix + ".fna")
    else:
        prefix = OUTPUT_DIR + "/{cls}/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}".format(cls=wildcards.cls,
         barcode=wildcards.barcode, c=wildcards.c, clust_id=wildcards.clust_id, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        OUTPUT_DIR + "/{cls}/{barcode}/{c}/split/{clust_id}.fastq",
        lambda wc: get_racon_input(wc),
    output: OUTPUT_DIR + "/{cls}/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}.fna"
    message: "Polish draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}] with racon, round={wildcards.iter} [{wildcards.cls}]"
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR +"/logs/{cls}/{barcode}/{c}/{clust_id}/racon/round{iter}.log"
    benchmark: OUTPUT_DIR +"/benchmarks/{cls}/{barcode}/{c}/{clust_id}/racon/round{iter}.txt"
    threads: config["threads"]["normal"]
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

rule medaka_consensus:
    input:
        fna = expand(OUTPUT_DIR + "/{{cls}}/{{barcode}}/{{c}}/polish/{{clust_id}}/draft/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = OUTPUT_DIR + "/{cls}/{barcode}/{c}/split/{clust_id}.fastq",
    output: 
        fasta = OUTPUT_DIR + "/{cls}/{barcode}/{c}/polish/{clust_id}/medaka/consensus.fasta",
        _dir = directory(OUTPUT_DIR + "/{cls}/{barcode}/{c}/polish/{clust_id}/medaka"),
    message: "Generate consensus in draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}] with medaka [{wildcards.cls}]"
    params:
        m = config["medaka"]["m"],
    conda: "../envs/polish.yaml"
    log: OUTPUT_DIR + "/logs/{cls}/{barcode}/{c}/{clust_id}/medaka.log"
    benchmark: OUTPUT_DIR + "/benchmarks/{cls}/{barcode}/{c}/{clust_id}/medaka.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        medaka_consensus -i {input.fastq} \
        -d {input.fna} -o {output._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        """

# isONcorCon
# pseudo checkpoint
checkpoint cls_isONcor:
    input: OUTPUT_DIR + "/isONclustCon/clusters"
    output: directory(OUTPUT_DIR + "/isONcorCon/clusters")
    shell: "cp -rf {input} {output}"

# run_isoncorrect provide multithread processing for isONcorrect in batches
# racon seems not to be supported in batch mode
rule isONcorrect:
    input: rules.get_fqs_split3.output,
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
 
# get polished asembly
# kmerCon
def get_kmerCon(wildcards):
    bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
    fnas = []
    for i in bc_kbs:
        bc, kb, = i.split("_")
        fnas.append(OUTPUT_DIR + "/kmerCon/{bc}/{kb}/polish/0/medaka/consensus.fasta".format(bc=bc, kb=kb))
    return fnas

rule collect_kmerCon:
    input: lambda wc: get_kmerCon(wc),
    output: OUTPUT_DIR + "/kmerCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i = [ i.split("/")[index] for index in [-6, -5] ]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "\n"
                        out.write(line)

# clustCon
def get_clustCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.cls_clustCon.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    fnas = []
    for i in bc_kb_cis:
        bc, kb, ci = i.split("_")
        fnas.append(OUTPUT_DIR + "/clustCon/{bc}/{kb}/polish/{ci}/medaka/consensus.fasta".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_clustCon:
    input: lambda wc: get_clustCon(wc),
    output: OUTPUT_DIR + "/clustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                barcode_i, c_i, id_i = [ i.split("/")[index] for index in [-6, -5, -3] ]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "\n"
                        out.write(line)

# isONclustCon
def get_isONclustCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.cls_isONclust.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    fnas = []
    for i in bc_kb_cis:
        bc, kb, ci = i.split("_")
        fnas.append(OUTPUT_DIR + "/isONclustCon/{bc}/{kb}/polish/{ci}/medaka/consensus.fasta".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_isONclustCon:
    input: lambda wc: get_isONclustCon(wc),
    output: OUTPUT_DIR + "/isONclustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                bc_i, kb_i, ci_i = [ i.split("/")[index] for index in [-6, -5, -3] ]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + bc_i + "_" + kb_i + "_" + ci_i + "\n"
                        out.write(line)

# isONcorCon   
def get_isONcorCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.cls_isONcor.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
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