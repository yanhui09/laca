def get_fq4Con(kmerbin = config["kmerbin"]):
    check_val("kmerbin", kmerbin, bool)
    if kmerbin == True:
        out = rules.split_bin.output
    else:
        out = rules.skip_bin.output
    return out

# kmerCon
rule get_fqs_split:
    input: get_fq4Con()
    output: temp("kmerCon/split/{barcode}_{c}_0.fastq"),
    log: "logs/kmerCon/{barcode}_{c}_0/get_fqs_split.log"
    benchmark: "benchmarks/kmerCon/{barcode}_{c}_0/get_fqs_split.txt"
    shell: "cp -f {input} {output} 2> {log}"

rule spoa:
    input: rules.get_fqs_split.output
    output: temp("kmerCon/polish/{barcode}_{c}_0/minimap2/raw.fna")
    conda: '../envs/spoa.yaml'
    params:
        l = config["spoa"]["l"],
        r = config["spoa"]["r"],
        g = config["spoa"]["g"],
    log: "logs/kmerCon/{barcode}_{c}_0/spoa.log"
    benchmark: "benchmarks/kmerCon/{barcode}_{c}_0/spoa.txt"
    shell: "spoa {input} -l {params.l} -r {params.r} -g {params.g} -s > {output} 2> {log}"

# clustCon
# draw draft with max average score from pairwise alignments
rule minimap2ava:
    input: get_fq4Con()
    output: temp("clustCon/minimap2va/{barcode}_{c}.paf")
    conda: '../envs/minimap2.yaml'
    params:
        x = config["minimap2"]["x_ava"],
    log: "logs/clustCon/minimap2ava/{barcode}_{c}.log"
    benchmark: "benchmarks/clustCon/minimap2ava/{barcode}_{c}.txt"
    threads: config["threads"]["large"]
    shell:
        "minimap2 -t {threads} -x {params.x} --no-long-join -r100"
        " {input} {input} > {output} 2> {log}"

rule ava2clust:
    input: rules.minimap2ava.output
    output: temp("clustCon/ava2clust/{barcode}_{c}.csv")
    conda: "../envs/clustCon.yaml"
    params:
        prefix = "clustCon/ava2clust/{barcode}_{c}",
        min_score_frac = config["ava2clust"]["min_score_frac"],
        min_reads = config["min_cluster_size"],
        max_recurs = config["ava2clust"]["max_recursion"],
    log: "logs/clustCon/ava2clust/{barcode}_{c}.log"
    benchmark: "benchmarks/clustCon/ava2clust/{barcode}_{c}.txt"
    shell:
        "python {workflow.basedir}/scripts/binClust.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input} > {log} 2>& 1"

def get_clust(wildcards, pool = config["pool"], kmerbin = config["kmerbin"]):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
        
    if kmerbin == True:
        bin2cls = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            bin2cls.append("clustCon/ava2clust/{bc}_{kb}.csv".format(bc=bc, kb=kb))
    else:
        if pool == True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        bin2cls = expand("clustCon/ava2clust/{bc}_all.csv", bc=bcs)
    return bin2cls

checkpoint cls_clustCon:
    input:
        ["kmerBin/clusters", ".qc_DONE"] if config["kmerbin"] else ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demultiplexed(wc)),
        lambda wc: get_kmerBin(wc),
        cls = lambda wc: get_clust(wc),
    output: directory("clustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.cls):
            num_lines = sum(1 for line in open(i))
            if num_lines > params.min_size:
                barcode, c = [ i.split("/")[-1].removesuffix(".csv").split("_")[index] for index in [-2, -1] ]
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
        binned = get_fq4Con(),
    output:
        pool = temp("clustCon/split/{barcode}_{c}_{clust_id}.fastq"),
        ref = temp("clustCon/polish/{barcode}_{c}_{clust_id}/minimap2/raw.fna"),
    conda: '../envs/seqkit.yaml'
    log: "logs/clustCon/{barcode}_{c}_{clust_id}/get_fqs_split.log"
    benchmark: "benchmarks/clustCon/{barcode}_{c}_{clust_id}/get_fqs_split.txt"
    shell:
        """
        seqkit grep -f {input.pool} {input.binned} -o {output.pool} --quiet 2> {log}
        seqkit grep -f {input.ref} {input.binned} --quiet | seqkit fq2fa | \
        seqkit replace -p '^(.+)$' -r '{wildcards.barcode}_{wildcards.c}_{wildcards.clust_id}' -o {output.ref} 2> {log}
        """

# isONclust
rule isONclust:
    input: get_fq4Con()
    output:
        _dir = temp(directory("isONclustCon/isONclust/{barcode}_{c}")),
        tsv = temp("isONclustCon/isONclust/{barcode}_{c}.tsv"),
    params:
        k = config["isONclust"]["k"],
        w = config["isONclust"]["w"],
    conda: "../envs/isONcorCon.yaml"
    log: "logs/isONclustCon/isONclust/{barcode}_{c}.log"
    benchmark: "benchmarks/isONclustCon/isONclust/{barcode}_{c}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        isONclust --k {params.k} --w {params.w} --fastq {input} --outfolder {output._dir} --t {threads} > {log} 2>&1
        mv {output._dir}/final_clusters.tsv {output.tsv}
        """

def get_isONclust(wildcards, pool = config["pool"], kmerbin = config["kmerbin"]):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
        
    if kmerbin == True:
        bin2cls = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            bin2cls.append("isONclustCon/isONclust/{bc}_{kb}.tsv".format(bc=bc, kb=kb))
    else:
        if pool == True:
           bcs = ["pooled"]
        else:
           bcs = get_qced(wildcards)
        bin2cls = expand("isONclustCon/isONclust/{bc}_all.tsv", bc=bcs)
    return bin2cls

checkpoint cls_isONclust:
    input:
        ["kmerBin/clusters", ".qc_DONE"] if config["kmerbin"] else ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demultiplexed(wc)),
        lambda wc: get_kmerBin(wc),
        cls = lambda wc: get_isONclust(wc),
    output: directory("isONclustCon/clusters")
    params:
        min_size = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.cls):
            barcode, c = [ i.split('/')[-1].removesuffix(".tsv").split('_')[index] for index in [-2, -1] ]
            df_i = pd.read_csv(i, sep = '\t', header = None)
            df_i.columns = ['cluster', 'seqid']
            for clust_id, df_clust in df_i.groupby('cluster'):
                if len(df_clust) >= params.min_size:
                    df_clust['seqid'].to_csv(output[0] + "/{barcode}_{c}_{clust_id}.csv".format(barcode=barcode, c=c, clust_id=clust_id),
                     header = False, index = False)

rule get_fqs_split3:
    input:
        bin2clust = "isONclustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        binned = get_fq4Con(),
    output: temp("isONclustCon/split/{barcode}_{c}_{clust_id}.fastq"),
    conda: '../envs/seqkit.yaml'
    log: "logs/isONclustCon/{barcode}_{c}_{clust_id}/get_fqs_split.log"
    benchmark: "benchmarks/isONclustCon/{barcode}_{c}_{clust_id}/get_fqs_split.txt"
    shell: "seqkit grep -f {input.bin2clust} {input.binned} -o {output} --quiet 2> {log}"

use rule spoa as spoa2 with:
    input:
        rules.get_fqs_split3.output
    output: 
        temp("isONclustCon/polish/{barcode}_{c}_{clust_id}/minimap2/raw.fna")
    log: 
        "logs/isONclustCon/{barcode}_{c}_{clust_id}/spoa.log"
    benchmark: 
        "benchmarks/isONclustCon/{barcode}_{c}_{clust_id}/spoa.txt"

# polish with racon and medaka
# reused in racon iterations
rule minimap2polish:
    input: 
      ref = "{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/{assembly}.fna",
      fastq = "{cls}/split/{barcode}_{c}_{clust_id}.fastq",
    output: temp("{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/{assembly}.paf"),
    message: "Polish draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}]: alignments against {wildcards.assembly} assembly [{wildcards.cls}]"
    params:
        x = config["minimap2"]["x_map"]
    conda: "../envs/minimap2.yaml"
    log: "logs/{cls}/{barcode}_{c}_{clust_id}/minimap2_{assembly}.log"
    benchmark: "benchmarks/{cls}/{barcode}_{c}_{clust_id}/minimap2_{assembly}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        # if ref is empty, make dummy output
        if [ ! -s {input.ref} ]; then
            touch {output} 2> {log}
        else
            minimap2 -t {threads} -x {params.x} {input.ref} {input.fastq} > {output} 2> {log}
        fi
        """

def get_racon_input(wildcards):
    # adjust input based on racon iteritions
    if int(wildcards.iter) == 1:
        prefix = "{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/raw"
    else:
        prefix = "{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/racon_{iter}".format(cls=wildcards.cls,
         barcode=wildcards.barcode, c=wildcards.c, clust_id=wildcards.clust_id, iter=str(int(wildcards.iter) - 1))
    return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        "{cls}/split/{barcode}_{c}_{clust_id}.fastq",
        lambda wc: get_racon_input(wc),
    output: temp("{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/racon_{iter}.fna")
    message: "Polish draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}] with racon, round={wildcards.iter} [{wildcards.cls}]"
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/racon.yaml"
    log: "logs/{cls}/{barcode}_{c}_{clust_id}/racon_{iter}.log"
    benchmark: "benchmarks/{cls}/{barcode}_{c}_{clust_id}/racon_{iter}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        # if paf file is empty, make dummy output
        if [ ! -s {input[1]} ]; then
            touch {output} 2> {log}
        else
            racon -m {params.m} -x {params.x} -g {params.g} -w {params.w} -t {threads} {input} > {output} 2> {log}
        fi
        """

# add iter for medaka
def get_medaka_files(wildcards, racon_iter = config["racon"]["iter"], index = False):
    if int(wildcards.iter2) == 1:
        fna = "{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/racon_{iter}.fna".format(
            cls=wildcards.cls, barcode=wildcards.barcode, c=wildcards.c, clust_id=wildcards.clust_id, iter=racon_iter)
    else:
        fna = "{cls}/polish/{barcode}_{c}_{clust_id}/medaka_{iter}/consensus.fasta".format(
            cls=wildcards.cls, barcode=wildcards.barcode, c=wildcards.c, clust_id=wildcards.clust_id, iter=str(int(wildcards.iter2) - 1))
    if index == True:
        return(fna + ".fai", fna + ".map-ont.mmi")
    else:
        return(fna)

rule medaka_consensus:
    input:
        fna = lambda wc: get_medaka_files(wc),
        fastq = "{cls}/split/{barcode}_{c}_{clust_id}.fastq",
    output: 
        temp(expand("{{cls}}/polish/{{barcode}}_{{c}}_{{clust_id}}/medaka_{{iter2}}/consensus{ext}",
        ext = [".fasta", ".fasta.gaps_in_draft_coords.bed", "_probs.hdf"])),
        temp(expand("{{cls}}/polish/{{barcode}}_{{c}}_{{clust_id}}/medaka_{{iter2}}/calls{ext}",
        ext = ["_to_draft.bam", "_to_draft.bam.bai"])),
    message: "Generate consensus in draft [barcode={wildcards.barcode}, bin={wildcards.c}, id={wildcards.clust_id}] with medaka, round={wildcards.iter2} [{wildcards.cls}]"
    params:
        m = config["medaka"]["m"],
        _dir = "{cls}/polish/{barcode}_{c}_{clust_id}/medaka_{iter2}",
        inedxs = lambda wc: get_medaka_files(wc, index = True),
    conda: "../envs/medaka.yaml"
    log: "logs/{cls}/{barcode}_{c}_{clust_id}/medaka_{iter2}.log"
    benchmark: "benchmarks/{cls}/{barcode}_{c}_{clust_id}/medaka_{iter2}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        # if fna file is empty, make dummy output
        if [ ! -s {input.fna} ]; then
            mkdir -p {params._dir} 2> {log}
            touch {output} 2>> {log}
        else
            export TF_FORCE_GPU_ALLOW_GROWTH=true
            medaka_consensus -i {input.fastq} -d {input.fna} -o {params._dir} -t {threads} -m {params.m} > {log} 2>&1
            rm -f {params.inedxs}
        fi
        """

# isONcorCon
# run_isoncorrect provide multithread processing for isONcorrect in batches
# racon seems not to be supported in batch mode
rule isONcorrect:
    input: rules.get_fqs_split3.output
    output: temp("isONcorCon/polish/{barcode}_{c}_{clust_id}/isONcor/{clust_id}/corrected_reads.fastq")
    conda: "../envs/isONcorCon.yaml"
    params:
        _dir = "isONcorCon/polish/{barcode}_{c}_{clust_id}",
        clust_id = "{clust_id}",
        max_seqs = 2000,
    log: "logs/isONcorCon/isONcorrect/{barcode}_{c}_{clust_id}.log"
    benchmark: "benchmarks/isONcorCon/isONcorrect/{barcode}_{c}_{clust_id}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        mkdir -p {params._dir}/isONclust
        cp {input} {params._dir}/isONclust/{params.clust_id}.fastq
        # number of reads in fastqs
        nlines=$(cat {params._dir}/isONclust/{params.clust_id}.fastq | wc -l)
        nreads=$((nlines / 4))
        if [ $nreads -gt {params.max_seqs} ]; then
            run_isoncorrect --t {threads} --fastq_folder {params._dir}/isONclust --outfolder {params._dir}/isONcor --set_w_dynamically --split_wrt_batches --max_seqs {params.max_seqs} > {log} 2>&1
        else
            run_isoncorrect --t {threads} --fastq_folder {params._dir}/isONclust --outfolder {params._dir}/isONcor --set_w_dynamically --max_seqs {params.max_seqs} > {log} 2>&1
        fi
        rm -rf {params._dir}/isONclust
        find {params._dir}/isONcor/{params.clust_id} -mindepth 1 ! -name 'corrected_reads.fastq' | xargs rm -rf
        """

rule isoCon:
    input: rules.isONcorrect.output
    output:
        cls = "isONcorCon/clusters/{barcode}_{c}_{clust_id}.tsv",
        fna = temp("isONcorCon/polish/{barcode}_{c}_{clust_id}/IsoCon/candidates.fasta"),
    conda: "../envs/isONcorCon.yaml"
    params:
        prefix = "isONcorCon/polish/{barcode}_{c}_{clust_id}",
    log: "logs/isONcorCon/isoCon/{barcode}_{c}_{clust_id}.log"
    benchmark: "benchmarks/isONcorCon/isoCon/{barcode}_{c}_{clust_id}.txt"
    threads: config["threads"]["large"]
    shell: 
        """
        IsoCon pipeline -fl_reads {input} -outfolder {params.prefix}/IsoCon --nr_cores {threads} > {log} 2>&1
        mv {params.prefix}/IsoCon/cluster_info.tsv {output.cls} 
        mv {params.prefix}/IsoCon/final_candidates.fa {output.fna}
        find {params.prefix}/IsoCon -mindepth 1 ! -name 'candidates.fasta' | xargs rm -rf 
        """
 
# get polished asembly
def merge_consensus(fi, fo):
    with open(fo, "w") as out:
        for i in fi:
            # if fi is empty, skip; rm dummy output
            if os.stat(i).st_size == 0:
                continue
            barcode_i, c_i, id_i = [ i.split("/")[-3].split("_")[index] for index in [-3, -2, -1] ]
            with open(i, "r") as inp:
                j = 1
                for line in inp:
                    if line.startswith(">"):
                        line = ">" + barcode_i + "_" + c_i + "_" + id_i + "_cand" + str(j) + "\n"
                        j += 1
                    out.write(line)

def get_clusters(wildcards):
    if wildcards.cls == "kmerCon":
        return "kmerBin/clusters"
    if wildcards.cls == "isONcorCon":
        return "isONclustCon/clusters"
    else:
        return "{cls}/clusters"

def get_consensus(wildcards, medaka_iter = config["medaka"]["iter"]):
    if wildcards.cls == "kmerCon":
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        bc_kb_cis = [i + "_0" for i in bc_kbs]
    elif wildcards.cls == "clustCon":
        bc_kb_cis = glob_wildcards(checkpoints.cls_clustCon.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    elif wildcards.cls == "isONclustCon" or wildcards.cls == "isONcorCon":
        bc_kb_cis = glob_wildcards(checkpoints.cls_isONclust.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    elif wildcards.cls == "umiCon":
        bc_kb_cis = glob_wildcards(checkpoints.cls_umiCon.get(**wildcards).output[0] + "/{bc_kb_ci}.txt").bc_kb_ci
    else:
        raise ValueError("Unknown consensus method: " + wildcards.cls)
    
    if wildcards.cls == "isONcorCon":
        return expand("{{cls}}/polish/{bc_kb_ci}/IsoCon/candidates.fasta", cls = wildcards.cls, bc_kb_ci = bc_kb_cis)
    else:
        return expand("{{cls}}/polish/{bc_kb_ci}/medaka_{iter}/consensus.fasta", cls = wildcards.cls, bc_kb_ci = bc_kb_cis, iter = medaka_iter)
            
rule collect_consensus:
    input: 
        lambda wc: get_clusters(wc),
        fna = lambda wc: get_consensus(wc),
    output: "{cls}/{cls}.fna"
    run: merge_consensus(fi = input.fna, fo = output[0]) 