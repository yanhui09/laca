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
    output: temp("kmerCon/{barcode}/{c}/split/0.fastq"),
    log: "logs/kmerCon/{barcode}/{c}/0/get_fqs_split.log"
    benchmark: "benchmarks/kmerCon/{barcode}/{c}/0/get_fqs_split.txt"
    shell: "cp -f {input} {output} 2> {log}"

rule spoa:
    input: rules.get_fqs_split.output
    output: temp("kmerCon/{barcode}/{c}/polish/0/draft/raw.fna")
    conda: '../envs/spoa.yaml'
    params:
        l = config["spoa"]["l"],
        r = config["spoa"]["r"],
        g = config["spoa"]["g"],
    log: "logs/isONclustCon/{barcode}/{c}/0/spoa_consensus.log"
    benchmark: "benchmarks/isONclustCon/{barcode}/{c}/0/spoa_consensus.txt"
    shell: "spoa {input} -l {params.l} -r {params.r} -g {params.g} -s > {output} 2> {log}"

# clustCon
# draw draft with max average score from pairwise alignments
rule minimap2clust:
    input: get_fq4Con(config["kmerbin"])
    output: temp("clustCon/{barcode}/avr_aln/{c}/minimap2clust.paf")
    conda: '../envs/minimap2.yaml'
    log: "logs/clustCon/{barcode}/minimap2clust/{c}.log"
    benchmark: "benchmarks/clustCon/{barcode}/minimap2clust/{c}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x ava-ont --no-long-join -r100"
        " {input} {input} > {output} 2> {log}"

rule bin2clustering:
    input: "clustCon/{barcode}/avr_aln/{c}/minimap2clust.paf"
    output: temp("clustCon/{barcode}/avr_aln/{c}/bin2clust.csv")
    conda: "../envs/clustCon.yaml"
    params:
        prefix = "clustCon/{barcode}/avr_aln/{c}/bin2clust",
        min_score_frac = config["clustCon"]["min_score_frac"],
        min_reads = config["min_cluster_size"],
        max_recurs = config["clustCon"]["max_recursion"],
    log: "logs/clustCon/{barcode}/bin2clust/{c}.log"
    benchmark: "benchmarks/clustCon/{barcode}/bin2clust/{c}.txt"
    shell:
        "python {workflow.basedir}/scripts/binClust.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input} > {log} 2>& 1"

def get_clust(wildcards, pool = True, kmerbin = True):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
        
    if kmerbin == True:
        bin2cls = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            bin2cls.append("clustCon/{bc}/avr_aln/{kb}/bin2clust.csv".format(bc=bc, kb=kb))
    else:
        if pool == True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        bin2cls = expand("clustCon/{bc}/avr_aln/all/bin2clust.csv", bc=bcs)
    return bin2cls

checkpoint cls_clustCon:
    input:
        ["kmerBin/clusters","qc/qfilt/empty"] if config["kmerbin"] else "qc/qfilt/empty",
        lambda wc: get_kmerBin(wc, pool = config["pool"], kmerbin = config["kmerbin"]),
        cls = lambda wc: get_clust(wc, pool = config["pool"], kmerbin = config["kmerbin"]),
    output: directory("clustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import shutil
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.cls):
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
        pool = "clustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        ref = "clustCon/clusters/{barcode}_{c}_{clust_id}.ref",
        binned = get_fq4Con(config["kmerbin"]),
    output:
        pool = temp("clustCon/{barcode}/{c}/split/{clust_id}.fastq"),
        ref = temp("clustCon/{barcode}/{c}/polish/{clust_id}/draft/raw.fna"),
    conda: '../envs/seqkit.yaml'
    log: "logs/clustCon/{barcode}/{c}/{clust_id}/get_fqs_split.log"
    benchmark: "benchmarks/clustCon/{barcode}/{c}/{clust_id}/get_fqs_split.txt"
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
        _dir = directory("isONclustCon/{barcode}/{c}"),
        tsv = "isONclustCon/{barcode}/{c}/final_clusters.tsv"
    conda: "../envs/isONcorCon.yaml"
    log: "logs/isONclustCon/{barcode}/{c}/isONclust.log"
    benchmark: "benchmarks/isONclustCon/{barcode}/{c}/isONclust.txt"
    threads: config["threads"]["large"]
    shell:
        "isONclust --ont --fastq {input} "
        "--outfolder {output._dir} --t {threads} > {log} 2>&1"

def get_isONclust(wildcards, pool = True, kmerbin = True):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
        
    if kmerbin == True:
        bin2cls = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            bin2cls.append("isONclustCon/{bc}/{kb}/final_clusters.tsv".format(bc=bc, kb=kb))
    else:
        if pool == True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        bin2cls = expand("isONclustCon/{bc}/all/final_clusters.tsv", bc=bcs)
    return bin2cls

checkpoint cls_isONclust:
    input:
        ["kmerBin/clusters","qc/qfilt/empty"] if config["kmerbin"] else "qc/qfilt/empty",
        lambda wc: get_kmerBin(wc, pool = config["pool"], kmerbin = config["kmerbin"]),
        cls = lambda wc: get_isONclust(wc, pool = config["pool"], kmerbin = config["kmerbin"]),
    output: directory("isONclustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.cls):
            barcode, c = [ i.split('/')[index] for index in [-3, -2] ]
            df_i = pd.read_csv(i, sep = '\t', header = None)
            df_i.columns = ['cluster', 'seqid']
            for clust_id, df_clust in df_i.groupby('cluster'):
                if len(df_clust) >= params.min_size:
                    df_clust['seqid'].to_csv(output[0] + "/{barcode}_{c}_{clust_id}.csv".format(barcode=barcode, c=c, clust_id=clust_id),
                     header = False, index = False)

rule get_fqs_split3:
    input:
        bin2clust = "isONclustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        binned = get_fq4Con(config["kmerbin"]),
    output: temp("isONclustCon/{barcode}/{c}/split/{clust_id}.fastq"),
    conda: '../envs/seqkit.yaml'
    log: "logs/isONclustCon/{barcode}/{c}/{clust_id}/get_fqs_split.log"
    benchmark: "benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/get_fqs_split.txt"
    shell: "seqkit grep -f {input.bin2clust} {input.binned} -o {output} --quiet 2> {log}"

use rule spoa as spoa2 with:
    input:
        rules.get_fqs_split3.output
    output: 
        temp("isONclustCon/{barcode}/{c}/polish/{clust_id}/draft/raw.fna")
    log: 
        "logs/isONclustCon/{barcode}/{c}/{clust_id}/spoa.log"
    benchmark: 
        "benchmarks/isONclustCon/{barcode}/{c}/{clust_id}/spoa.txt"

# polish with racon and medaka
# reused in racon iterations
rule minimap2polish:
    input: 
      ref = "{cls}/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.fna",
      fastq = "{cls}/{barcode}/{c}/split/{clust_id}.fastq",
    output: temp("{cls}/{barcode}/{c}/polish/{clust_id}/draft/{assembly}.paf"),
    message: "Polish draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}]: alignments against {wildcards.assembly} assembly [{wildcards.cls}]"
    params:
        x = config["minimap"]["x"]
    conda: "../envs/minimap2.yaml"
    log: "logs/{cls}/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.log"
    benchmark: "benchmarks/{cls}/{barcode}/{c}/{clust_id}/minimap2polish/{assembly}.txt"
    threads: config["threads"]["normal"]
    shell:
        "minimap2 -t {threads} -x {params.x}"
        " {input.ref} {input.fastq} > {output} 2> {log}"

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = "{cls}/{barcode}/{c}/polish/{clust_id}/draft/raw"
        return(prefix + ".paf", prefix + ".fna")
    else:
        prefix = "{cls}/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}".format(cls=wildcards.cls,
         barcode=wildcards.barcode, c=wildcards.c, clust_id=wildcards.clust_id, iter=str(int(wildcards.iter) - 1))
        return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        "{cls}/{barcode}/{c}/split/{clust_id}.fastq",
        lambda wc: get_racon_input(wc),
    output: temp("{cls}/{barcode}/{c}/polish/{clust_id}/draft/racon_{iter}.fna")
    message: "Polish draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}] with racon, round={wildcards.iter} [{wildcards.cls}]"
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/racon.yaml"
    log: "logs/{cls}/{barcode}/{c}/{clust_id}/racon/round{iter}.log"
    benchmark: "benchmarks/{cls}/{barcode}/{c}/{clust_id}/racon/round{iter}.txt"
    threads: config["threads"]["normal"]
    shell:
        "racon -m {params.m} -x {params.x}"
        " -g {params.g} -w {params.w} -t {threads}"
        " {input} > {output} 2> {log}"

rule medaka_consensus:
    input:
        fna = expand("{{cls}}/{{barcode}}/{{c}}/polish/{{clust_id}}/draft/racon_{iter}.fna", 
        iter = config["racon"]["iter"]),
        fastq = "{cls}/{barcode}/{c}/split/{clust_id}.fastq",
    output: 
        temp(expand("{{cls}}/{{barcode}}/{{c}}/polish/{{clust_id}}/medaka/consensus{ext}",
        ext = [".fasta", ".fasta.gaps_in_draft_coords.bed", "_probs.hdf"])),
        temp(expand("{{cls}}/{{barcode}}/{{c}}/polish/{{clust_id}}/medaka/calls{ext}",
        ext = ["_to_draft.bam", "_to_draft.bam.bai"])),
        temp(expand("{{cls}}/{{barcode}}/{{c}}/polish/{{clust_id}}/draft/racon_{iter}.fna{ext}", 
        iter = config["racon"]["iter"], ext = [".fai", ".map-ont.mmi"])),
    message: "Generate consensus in draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}] with medaka [{wildcards.cls}]"
    params:
        m = config["medaka"]["m"],
        _dir = "{cls}/{barcode}/{c}/polish/{clust_id}/medaka",
    conda: "../envs/medaka.yaml"
    log: "logs/{cls}/{barcode}/{c}/{clust_id}/medaka.log"
    benchmark: "benchmarks/{cls}/{barcode}/{c}/{clust_id}/medaka.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        export TF_FORCE_GPU_ALLOW_GROWTH=true
        medaka_consensus -i {input.fastq} \
        -d {input.fna} -o {params._dir} \
        -t {threads} -m {params.m} > {log} 2>&1;
        """

# isONcorCon
# run_isoncorrect provide multithread processing for isONcorrect in batches
# racon seems not to be supported in batch mode
rule isONcorrect:
    input: rules.get_fqs_split3.output,
    output:
        _dir = directory("isONcorCon/{barcode}/{c}/{clust_id}/isONcor"), 
        fq = temp("isONcorCon/{barcode}/{c}/{clust_id}/isONcor/{clust_id}/corrected_reads.fastq")
    conda: "../envs/isONcorCon.yaml"
    params:
        _dir = "isONcorCon/{barcode}/{c}/{clust_id}",
        max_seqs = 2000,
    log: "logs/isONcorCon/{barcode}/{c}/{clust_id}/isONcorrect.log"
    benchmark: "benchmarks/isONcorCon/{barcode}/{c}/{clust_id}/isONcorrect.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        mkdir -p {params._dir}/isONclust
        cp {input} {params._dir}/isONclust
        # number of reads in fastqs
        nlines=$(cat {params._dir}/isONclust/* | wc -l)
        nreads=$((nlines / 4))
        if [ $nreads -gt {params.max_seqs} ]; then
            run_isoncorrect --t {threads} --fastq_folder {params._dir}/isONclust --outfolder {output._dir} \
            --set_w_dynamically --split_wrt_batches --max_seqs {params.max_seqs} > {log} 2>&1
        else
            run_isoncorrect --t {threads} --fastq_folder {params._dir}/isONclust --outfolder {output._dir} \
            --set_w_dynamically --max_seqs {params.max_seqs} > {log} 2>&1
        fi
        rm -rf {params._dir}/isONclust
        """

rule isoCon:
    input: rules.isONcorrect.output.fq
    output:
        _dir = directory("isONcorCon/{barcode}/{c}/{clust_id}/isoCon"),
        fna = temp("isONcorCon/{barcode}/{c}/{clust_id}/isoCon/final_candidates.fa"),
    conda: "../envs/isONcorCon.yaml"
    log: "logs/isONcorCon/{barcode}/{c}/{clust_id}/isoCon.log"
    benchmark: "benchmarks/isONcorCon/{barcode}/{c}/{clust_id}/isoCon.txt"
    threads: config["threads"]["normal"]
    shell: "IsoCon pipeline -fl_reads {input} -outfolder {output._dir} --nr_cores {threads} > {log} 2>&1"
 
# get polished asembly
# kmerCon
def get_kmerCon(wildcards):
    bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
    fnas = []
    for i in bc_kbs:
        bc, kb, = i.split("_")
        fnas.append("kmerCon/{bc}/{kb}/polish/0/medaka/consensus.fasta".format(bc=bc, kb=kb))
    return fnas

rule collect_kmerCon:
    input: 
        "kmerBin/clusters",
        fna = lambda wc: get_kmerCon(wc),
    output: "kmerCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input.fna:
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
        fnas.append("clustCon/{bc}/{kb}/polish/{ci}/medaka/consensus.fasta".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_clustCon:
    input:
        "clustCon/clusters",
        fna = lambda wc: get_clustCon(wc),
    output: "clustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input.fna:
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
        fnas.append("isONclustCon/{bc}/{kb}/polish/{ci}/medaka/consensus.fasta".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_isONclustCon:
    input:
        "isONclustCon/clusters",
        fna = lambda wc: get_isONclustCon(wc),
    output: "isONclustCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input.fna:
                bc_i, kb_i, ci_i = [ i.split("/")[index] for index in [-6, -5, -3] ]
                with open(i, "r") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + bc_i + "_" + kb_i + "_" + ci_i + "\n"
                        out.write(line)

# isONcorCon   
def get_isONcorCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.cls_isONclust.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    fnas = []
    for i in bc_kb_cis:
        bc, kb, ci = i.split("_")
        fnas.append("isONcorCon/{bc}/{kb}/{ci}/isoCon/final_candidates.fa".format(bc=bc, kb=kb, ci=ci))
    return fnas

rule collect_isONcorCon:
    input: 
        "isONclustCon/clusters",
        fna = lambda wc: get_isONcorCon(wc),
    output: "isONcorCon.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input.fna:
                barcode_i, c_i, id_i = [ i.split("/")[index] for index in [-5, -4, -3] ]
                with open(i, "r") as inp:
                    j = 0
                    for line in inp:
                        if line.startswith(">"):
                            line = ">" + barcode_i + "_" + c_i + "_" + id_i + "_cand" + str(j) + "\n"
                            j += 1
                        out.write(line)