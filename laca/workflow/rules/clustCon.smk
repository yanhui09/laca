localrules: fqs_split1, cls_clustCon, fqs_split2, cls_isONclustCon, fqs_split3, cls_isONcorCon, get_IsoCon_cand, fqs_split4, cls_isONclustCon2, fqs_split_isONclust, collect_consensus 
# kmerCon
use rule fqs_split as fqs_split1 with:
    input:
        cluster = "kmerBin/clusters/{barcode}_{c}.csv",
        fqs = "qc/qfilt/{barcode}.fastq",
    output: 
        temp("kmerCon/split/{barcode}_{c}_0.fastq"),
    log: 
        "logs/kmerCon/fqs_split/{barcode}_{c}.log"
    benchmark: 
        "benchmarks/kmerCon/fqs_split/{barcode}_{c}.txt"

# clustCon
# draw draft with max average score from pairwise alignments
rule minimap2ava:
    input: rules.fqs_split.output
    output: temp("clustCon/minimap2va/{barcode}_{c}.paf")
    conda: '../envs/minimap2.yaml'
    params:
        x = config["minimap2"]["x_ava"],
    log: "logs/clustCon/minimap2ava/{barcode}_{c}.log"
    benchmark: "benchmarks/clustCon/minimap2ava/{barcode}_{c}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
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
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell:
        "python {workflow.basedir}/scripts/binClust.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input} > {log} 2>& 1"

def get_clust(wildcards):
    bin2cls = []
    bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
    for i in bc_kbs:
        bc, kb = i.split("_")
        bin2cls.append("clustCon/ava2clust/{bc}_{kb}.csv".format(bc=bc, kb=kb))
    return bin2cls

checkpoint cls_clustCon:
    input:
        ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: get_kmerBin(wc, kmerbin=True),
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
                # cluster starts from 0
                df["cluster"] = df["cluster"] - 1
                for clust_id, df_clust in df.groupby('cluster'):
                    if len(df_clust) >= params.min_size:
                        bc_kb_ci = "/{barcode}_{c}_{clust_id}".format(barcode=barcode, c=c, clust_id=clust_id)
                        df_clust['read_id'].to_csv(output[0] + bc_kb_ci + ".csv", header = False, index = False)
                        ref_idx  = df_clust['clust_read_score'].idxmax()
                        ref_read = df_clust.loc[ref_idx, ['read_id']]
                        ref_read.to_csv(output[0] + bc_kb_ci + ".ref", header=False, index=False)
   
rule fqs_split2:
    input:
        pool = "clustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        ref = "clustCon/clusters/{barcode}_{c}_{clust_id}.ref",
        binned = rules.fqs_split.output,
    output:
        pool = temp("clustCon/split/{barcode}_{c}_{clust_id}.fastq"),
        ref = temp("clustCon/polish/{barcode}_{c}_{clust_id}/minimap2/raw.fna"),
    log: "logs/clustCon/{barcode}_{c}_{clust_id}/fqs_split.log"
    benchmark: "benchmarks/clustCon/{barcode}_{c}_{clust_id}/fqs_split.txt"
    shell:
        """
        seqkit grep -f {input.pool} {input.binned} -o {output.pool} --quiet 2> {log}
        seqkit grep -f {input.ref} {input.binned} --quiet | seqkit fq2fa | \
        seqkit replace -p '^(.+)$' -r '{wildcards.barcode}_{wildcards.c}_{wildcards.clust_id}' -o {output.ref} 2> {log}
        """

def get_fq4Con(kmerbin = config["kmerbin"]):
    check_val("kmerbin", kmerbin, bool)
    if kmerbin == True:
        out = rules.fqs_split.output
    else:
        out = "qc/qfilt/{barcode}.fastq"
    return out

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
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell:
        """
        isONclust --k {params.k} --w {params.w} --fastq {input} --outfolder {output._dir} --t {threads} > {log} 2>&1
        mv {output._dir}/final_clusters.tsv {output.tsv}
        """

def get_isONclust(wildcards, cls, pool = config["pool"], kmerbin = config["kmerbin"]):
    check_val("pool", pool, bool)
    check_val("kmerbin", kmerbin, bool)
        
    if kmerbin == True:
        bin2cls = []
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        for i in bc_kbs:
            bc, kb = i.split("_")
            bin2cls.append("{cls}/isONclust/{bc}_{kb}.tsv".format(cls=cls, bc=bc, kb=kb))
    else:
        if pool == True:
           bcs = ["pooled"]
        else:
           bcs = get_qced_barcodes(wildcards)
        bin2cls = expand(cls + "/isONclust/{bc}_all.tsv", bc=bcs)
    return bin2cls

checkpoint cls_isONclustCon:
    input:
        ["kmerBin/clusters", ".qc_DONE"] if config["kmerbin"] else ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: get_kmerBin(wc),
        cls = lambda wc: get_isONclust(wc, cls="isONclustCon"),
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

use rule fqs_split as fqs_split3 with:
    input: 
        cluster = "isONclustCon/clusters/{barcode}_{c}_{clust_id}.csv",
        fqs = get_fq4Con()
    output: 
        temp("isONclustCon/split/{barcode}_{c}_{clust_id}.fastq")
    log: 
        "logs/isONclustCon/fqs_split/{barcode}_{c}_{clust_id}.log"
    benchmark:
        "benchmarks/isONclustCon/fqs_split/{barcode}_{c}_{clust_id}.txt"

# isONcorCon
# run_isoncorrect provide multithread processing for isONcorrect in batches
rule isONcorrect:
    input: rules.fqs_split3.output
    output: temp("isONcorCon/{barcode}_{c}_{clust_id}/isONcor/{clust_id}/corrected_reads.fastq")
    conda: "../envs/isONcorCon.yaml"
    params:
        _dir = "isONcorCon/{barcode}_{c}_{clust_id}",
        clust_id = "{clust_id}",
        max_seqs = 2000 * 4,
    log: "logs/isONcorCon/isONcorrect/{barcode}_{c}_{clust_id}.log"
    benchmark: "benchmarks/isONcorCon/isONcorrect/{barcode}_{c}_{clust_id}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    shell:
        """
        mkdir -p {params._dir}/isONclust
        cp {input} {params._dir}/isONclust/{params.clust_id}.fastq
        # number of reads in fastqs
        nlines=$(cat {params._dir}/isONclust/{params.clust_id}.fastq | wc -l)
        if [ $nlines -gt {params.max_seqs} ]; then
            run_isoncorrect --t {threads} --fastq_folder {params._dir}/isONclust --outfolder {params._dir}/isONcor --set_w_dynamically --split_wrt_batches --max_seqs {params.max_seqs} > {log} 2>&1
        else
            run_isoncorrect --t {threads} --fastq_folder {params._dir}/isONclust --outfolder {params._dir}/isONcor --set_w_dynamically --max_seqs {params.max_seqs} > {log} 2>&1
        fi
        rm -rf {params._dir}/isONclust
        find {params._dir}/isONcor/{params.clust_id} -mindepth 1 ! -name 'corrected_reads.fastq' | xargs rm -rf
        """

def get_isoCon_input(isONcor = config["isONcor"]):
    check_val("isONcor", isONcor, bool)
    if isONcor == True:
        return rules.isONcorrect.output
    else:
        return rules.fqs_split3.output

def check_isoCon_batch(batch_size = config["IsoCon"]["max_batch_size"]):
    check_val("IsoCon batch_size", batch_size, int)
    if int(batch_size) < -1:
        raise ValueError("IsoCon batch_size only accepts integer >= -1.")
check_isoCon_batch()

rule isoCon:
    input: get_isoCon_input()
    output:
        cls = temp("isONcorCon/{barcode}_{c}_{clust_id}/IsoCon/cluster_info.tsv"),
        fna = temp("isONcorCon/{barcode}_{c}_{clust_id}/IsoCon/final_candidates.fa"),
    conda: "../envs/isONcorCon.yaml"
    params:
        prefix = "isONcorCon/{barcode}_{c}_{clust_id}",
        min_candidates = config["min_cluster_size"],
        neighbor_search_depth =  int(config["IsoCon"]["neighbor_search_depth"]) if config["IsoCon"]["neighbor_search_depth"] else 2**32,
        p_value_threshold = config["IsoCon"]["p_value_threshold"],
        max_batch_size = -1 if int(config["IsoCon"]["max_batch_size"]) == -1 else int(config["IsoCon"]["max_batch_size"]) * 4,
    log: "logs/isONcorCon/isoCon/{barcode}_{c}_{clust_id}.log"
    benchmark: "benchmarks/isONcorCon/isoCon/{barcode}_{c}_{clust_id}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["long"],
    shell: 
        """
        nlines=$(cat {input} | wc -l)
        if [ {params.max_batch_size} -eq -1 ] || [ $nlines -le {params.max_batch_size} ]; then
          IsoCon pipeline -fl_reads {input} -outfolder {params.prefix}/IsoCon --nr_cores {threads} \
          --neighbor_search_depth {params.neighbor_search_depth} --p_value_threshold {params.p_value_threshold} \
          --prefilter_candidates --min_candidate_support {params.min_candidates} --cleanup >{log} 2>&1
          find {params.prefix}/IsoCon -mindepth 1 ! -name 'final_candidates.fa' ! -name 'cluster_info.tsv' | xargs rm -rf
        else
          # split input fastq into batches
          mkdir -p {params.prefix}/IsoCon/batches
          # determine the minimum partion size (number of batches), ceiling division
          min_part_size=$(((nlines + {params.max_batch_size} - 1) / {params.max_batch_size}))
          # determine the number of lines per batch, ceiling division, the nearest multiples of 4
          nlines_per_batch=$(((nlines + min_part_size * 4 - 1) / (min_part_size * 4) * 4)) 
          split -l $nlines_per_batch -a2 -d --additional-suffix='.fastq' {input} {params.prefix}/IsoCon/batches/b >{log} 2>&1
          for fq in {params.prefix}/IsoCon/batches/b*.fastq; do
            batch_id=$(basename $fq | cut -d'.' -f1)
            if [ -f {params.prefix}/IsoCon/batches/$batch_id/final_candidates.fa ] && [ -f {params.prefix}/IsoCon/batches/$batch_id/cluster_info.tsv ]; then
              continue
            fi
            IsoCon pipeline -fl_reads $fq -outfolder {params.prefix}/IsoCon/batches/$batch_id --nr_cores {threads} \
            --neighbor_search_depth {params.neighbor_search_depth} --p_value_threshold {params.p_value_threshold} \
            --prefilter_candidates --min_candidate_support {params.min_candidates} --cleanup >> {log} 2>&1
            find {params.prefix}/IsoCon/batches/$batch_id -mindepth 1 ! -name 'final_candidates.fa' ! -name 'cluster_info.tsv' | xargs rm -rf
            sed -i "s/_support/${{batch_id}}_support/" {params.prefix}/IsoCon/batches/$batch_id/cluster_info.tsv
            sed -i "s/_support/${{batch_id}}_support/" {params.prefix}/IsoCon/batches/$batch_id/final_candidates.fa
          done
          # merge the results
          cat {params.prefix}/IsoCon/batches/b*/cluster_info.tsv > {params.prefix}/IsoCon/cluster_info.tsv
          cat {params.prefix}/IsoCon/batches/b*/final_candidates.fa > {params.prefix}/IsoCon/final_candidates.fa
          rm -rf {params.prefix}/IsoCon/batches
        fi
        """

def get_isONcorCon(wildcards):
    bc_kb_cis = glob_wildcards(checkpoints.cls_isONclustCon.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    return {
    "cls": expand("isONcorCon/{bc_kb_ci}/IsoCon/cluster_info.tsv", bc_kb_ci = bc_kb_cis),
    "cands": expand("isONcorCon/{bc_kb_ci}/IsoCon/final_candidates.fa", bc_kb_ci = bc_kb_cis)
    }

checkpoint cls_isONcorCon:
    input: 
        ["kmerBin/clusters", ".qc_DONE"] if config["kmerbin"] else ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: expand("isONclustCon/split/{bc_kb_ci}.fastq", bc_kb_ci=glob_wildcards(checkpoints.cls_isONclustCon.get(**wc).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci),
        unpack(get_isONcorCon),
    output: directory("isONcorCon/clusters"),
    params:
        min_candidates = config["min_cluster_size"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.cls):
            # if empty, skip
            if os.stat(i).st_size == 0:
                continue
            bc_kb_ci = i.split('/')[-3]
            df_i = pd.read_csv(i, sep = '\t', header = None, usecols=range(2))
            df_i.columns = ['seqid', 'cluster']
            # only leave candidate id 'transcript_id_support_num'
            df_i['cluster'] = df_i['cluster'].apply(lambda x: x.split('_')[1])
            for clust_id, df_clust in df_i.groupby('cluster'):
                if len(df_clust) < params.min_candidates:
                    continue
                df_clust['seqid'].to_csv(output[0] + "/{bc_kb_ci}cand{clust_id}.csv".format(bc_kb_ci=bc_kb_ci, clust_id=clust_id),
                header = False, index = False)

rule get_IsoCon_cand:
    input:
        cls = "isONcorCon/clusters/{barcode}_{c}_{clust_id}cand{cand}.csv",
        cands = rules.isoCon.output.fna,
    output: temp("isONcorCon/polish/{barcode}_{c}_{clust_id}cand{cand}/minimap2/raw.fna")
    run:
            outdir = os.path.dirname(output[0])
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            bc_kb_ci, cand_id = input.cls.removesuffix('.csv').split('/')[-1].split('cand')
            with open (input.cands, 'r') as fi:
                lines = fi.readlines()

            for i, line in enumerate(lines):
                if line.startswith('>transcript_' + cand_id + '_support'):
                    header = '>' + bc_kb_ci + 'cand' + cand_id + '\n'
                    with open (output[0], 'w') as fo:
                        fo.write(header)
                        fo.write(lines[i+1])
                    break

rule fqs_split4:
    input:
        bin2clust = "isONcorCon/clusters/{barcode}_{c}_{clust_id}cand{cand}.csv",
        binned = rules.fqs_split3.output,
    output: temp("isONcorCon/split/{barcode}_{c}_{clust_id}cand{cand}.fastq")
    log: "logs/isONcorCon/{barcode}_{c}_{clust_id}cand{cand}/fqs_split.log"
    benchmark: "benchmarks/isONcorCon/{barcode}_{c}_{clust_id}cand{cand}/fqs_split.txt"
    shell: "seqkit grep -f {input.bin2clust} {input.binned} -o {output} --quiet 2> {log}"

use rule kmer_freqs as kmer_freqs2 with:
    input: 
        "{cls}/split/{barcode}_{c}_{clust_id}.fastq"
    output: 
        temp("{cls}/kmerBin/{barcode}_{c}_{clust_id}/kmer_freqs.txt")
    log: 
        "logs/{cls}/kmerBin/kmer_freqs/{barcode}_{c}_{clust_id}.log"
    benchmark: 
        "benchmarks/{cls}/kmerBin/kmer_freqs/{barcode}_{c}_{clust_id}.txt"

use rule umap as umap2 with:
    input: 
        rules.kmer_freqs2.output
    output: 
        cluster="{cls}/kmerBin/{barcode}_{c}_{clust_id}/hdbscan.tsv",
	    plot="{cls}/kmerBin/{barcode}_{c}_{clust_id}/hdbscan.png",
    log: 
        "logs/{cls}/kmerBin/umap/{barcode}_{c}_{clust_id}.log"
    benchmark: 
        "benchmarks/{cls}/kmerBin/umap/{barcode}_{c}_{clust_id}.txt"

checkpoint cls_isONclustCon2:
    input:
        ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: expand("isONclustCon/split/{bc_kb_ci}.fastq", bc_kb_ci=glob_wildcards(checkpoints.cls_isONclustCon.get(**wc).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci),
        bin = lambda wc: expand("isONclustCon/kmerBin/{bc_kb_ci}/hdbscan.tsv", bc_kb_ci=glob_wildcards(checkpoints.cls_isONclustCon.get(**wc).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci),
    output: directory("isONclustCon2/clusters"),
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.bin):
            bc_kb_ci = i.split("/")[-2]
            df_i = pd.read_csv(i, sep="\t")
            # unclustered reads are assigned to bin -1, drop
            df_i = df_i[df_i.bin_id != -1]

            if 'batch_id' in df_i.columns:
                 # concatanate the batch_id bin_id column
                 df_i['bin_id'] = df_i['bin_id'].astype(str) + df_i['batch_id'].astype(str)
            df_i = df_i[["read", "bin_id"]]
            for clust_id, df_clust in df_i.groupby('bin_id'):
                df_clust['read'].to_csv(output[0] + "/{bc_kb_ci}cand{c}.csv".format(
                    bc_kb_ci=bc_kb_ci, c=clust_id), header = False, index = False)

use rule fqs_split as fqs_split_isONclust with:
    input:
        cluster = "isONclustCon2/clusters/{barcode}_{c}_{clust_id}cand{cand}.csv",
        fqs = "isONclustCon/split/{barcode}_{c}_{clust_id}.fastq"
    output: 
        temp("isONclustCon2/split/{barcode}_{c}_{clust_id}cand{cand}.fastq"),
    log: 
        "logs/isONclustCon2/fqs_split/{barcode}_{c}_{clust_id}cand{cand}.log"
    benchmark: 
        "benchmarks/isONclustCon2/fqs_split/{barcode}_{c}_{clust_id}cand{cand}.txt"

rule spoa:
    input: "{cls}/split/{barcode}_{c}_{clust_id}.fastq"
    output: temp("{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/raw.fna")
    conda: '../envs/isONcorCon.yaml'
    params:
        l = config["spoa"]["l"],
        r = config["spoa"]["r"],
        g = config["spoa"]["g"],
        M = config["spoa"]["M"],
    log: "logs/{cls}/spoa/{barcode}_{c}_{clust_id}.log"
    benchmark: "benchmarks/{cls}/spoa/{barcode}_{c}_{clust_id}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["long"],
    shell: 
        """
        # touch if not empty
        if [ -s {output} ]; then
            touch {output}
        else
            # use first 4*M reads if M > 0
            if [ {params.M} -gt 0 ]; then
                head -n $((4*{params.M})) {input} > {input}.M.fastq
                spoa {input}.M.fastq -l {params.l} -r {params.r} -g {params.g} -s > {output} 2> {log}
                rm {input}.M.fastq
            else
                spoa {input} -l {params.l} -r {params.r} -g {params.g} -s > {output} 2> {log}
            fi    
        fi
        """

ruleorder: fqs_split2 > spoa
ruleorder: get_IsoCon_cand > spoa

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
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
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
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
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
        if int(racon_iter) == 0:
            fna = "{cls}/polish/{barcode}_{c}_{clust_id}/minimap2/raw.fna".format(
                cls=wildcards.cls, barcode=wildcards.barcode, c=wildcards.c, clust_id=wildcards.clust_id)
        else:
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
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
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

# get polished asembly
def merge_consensus(fi, fo):
    with open(fo, "w") as out:
        for i in fi:
            # if fi is empty, skip; rm dummy output
            if os.stat(i).st_size == 0:
                continue
            bc_kb_ci_cand = i.split("/")[-3]
            if "cand" not in bc_kb_ci_cand:
                bc_kb_ci_cand = bc_kb_ci_cand + "cand1"
            
            with open(i, "r") as inp:
                for line in inp:
                    if line.startswith(">"):
                        line = ">" + bc_kb_ci_cand + "\n"
                    out.write(line)

def get_clusters(wildcards):
    if wildcards.cls == "kmerCon":
        return "kmerBin/clusters"
    else:
        return "{cls}/clusters"

def get_consensus(wildcards, medaka_iter = config["medaka"]["iter"], racon_iter = config["racon"]["iter"]):
    # iter >= 0, integer
    check_val("racon iter", racon_iter, int)
    check_val("medaka iter", medaka_iter, int)
    if racon_iter < 0 or medaka_iter < 0:
        raise ValueError("racon and medaka iter shall be >= 0")

    if wildcards.cls == "kmerCon":
        bc_kbs = glob_wildcards(checkpoints.cls_kmerbin.get(**wildcards).output[0] + "/{bc_kb}.csv").bc_kb
        candidates = [i + "_0" for i in bc_kbs]
    elif wildcards.cls == "clustCon":
        candidates = glob_wildcards(checkpoints.cls_clustCon.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    elif wildcards.cls == "isONclustCon":
        candidates = glob_wildcards(checkpoints.cls_isONclustCon.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    elif wildcards.cls == "isONclustCon2":
        candidates = glob_wildcards(checkpoints.cls_isONclustCon2.get(**wildcards).output[0] + "/{bc_kb_ci_cand}.csv").bc_kb_ci_cand
    elif wildcards.cls == "isONcorCon":
        candidates = glob_wildcards(checkpoints.cls_isONcorCon.get(**wildcards).output[0] + "/{bc_kb_ci_cand}.csv").bc_kb_ci_cand
    elif wildcards.cls == "umiCon":
        candidates = glob_wildcards(checkpoints.cls_umiCon.get(**wildcards).output[0] + "/{bc_kb_ci}.csv").bc_kb_ci
    else:
        raise ValueError("Unknown consensus method: " + wildcards.cls)
    
    if medaka_iter == 0:
        if racon_iter == 0:
            return expand("{{cls}}/polish/{cand}/minimap2/raw.fna", cls = wildcards.cls, cand = candidates)
        else:
            return expand("{{cls}}/polish/{cand}/minimap2/racon_{iter}.fna", cls = wildcards.cls, cand = candidates, iter = racon_iter)
    else:
        return expand("{{cls}}/polish/{cand}/medaka_{iter}/consensus.fasta", cls = wildcards.cls, cand = candidates, iter = medaka_iter)
            
rule collect_consensus:
    input: 
        lambda wc: get_clusters(wc),
        fna = lambda wc: get_consensus(wc),
    output: "{cls}/{cls}.fna"
    run: merge_consensus(fi = input.fna, fo = output[0])