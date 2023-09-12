localrules: fqs_split_kmerCon, cls_miniCon, fqs_split_miniCon, cls_isoCon, cls_isoCon, get_IsoCon_cand, fqs_split_isoCon, collect_consensus 
# kmerCon
rule fqs_split_kmerCon:
    input:
        members = rules.fqs_split_meshclust.output.members,
        centroid = rules.fqs_split_meshclust.output.centroid,
    output: 
        members = temp("kmerCon/split/{barcode}_{c1}_{c2}_{c3}cand1.fastq"),
        centroid = temp("kmerCon/polish/{barcode}_{c1}_{c2}_{c3}cand1/minimap2/raw.fna"),
    shell: "cp {input.members} {output.members} && cp {input.centroid} {output.centroid}"

def get_kmerCon(wildcards):
    bc_clss = glob_wildcards(checkpoints.cls_meshclust.get(**wildcards).output[0] + "/{bc_cls}.csv").bc_cls
    return {
    "fastq": expand("kmerCon/split/{bc_cls}cand1.fastq", bc_cls = bc_clss),
    "fna": expand("kmerCon/polish/{bc_cls}cand1/minimap2/raw.fna", bc_cls = bc_clss),
    "clss": expand("clust/clusters/{bc_cls}.csv", bc_cls = bc_clss),
    }

checkpoint cls_kmerCon:
    input:
        "clust/clusters", ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: get_kmerclust(wc),
        unpack(get_kmerCon),
    output: directory("kmerCon/clusters")
    params:
        min_size = config["min_support_reads"],
    run:
        import shutil

        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.clss):
            num_lines = sum(1 for line in open(i))
            if num_lines > params.min_size:
                bc_clss = i.split("/")[-1].removesuffix(".csv")
                bc_clss_cand = "/{bc_clss}cand1".format(bc_clss=bc_clss)
                shutil.copyfile(i, output[0] + bc_clss_cand + ".csv")

# miniCon
# refine cluster with max average score from base-level alignments
rule minimap2ava:
    input: rules.fqs_split_meshclust.output.members
    output: temp("miniCon/minimap2ava/{barcode}_{c1}_{c2}_{c3}.paf")
    conda: "../envs/minimap2.yaml"
    params:
        x = config["minimap2"]["x_ava"],
        f = config["yacrd"]["minimap2"]["f"],
    log: "logs/miniCon/minimap2ava/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/miniCon/minimap2ava/{barcode}_{c1}_{c2}_{c3}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell:
        "minimap2 -t {threads} -x {params.x} -f {params.f} --no-long-join -r100"
        " {input} {input} > {output} 2> {log}"

rule ava2clust:
    input: rules.minimap2ava.output
    output: temp("miniCon/ava2clust/{barcode}_{c1}_{c2}_{c3}.csv")
    conda: "../envs/miniCon.yaml"
    params:
        prefix = "miniCon/ava2clust/{barcode}_{c1}_{c2}_{c3}",
        min_score_frac = config["ava2clust"]["min_score_frac"],
        min_reads = config["min_support_reads"],
        max_recurs = config["ava2clust"]["max_recursion"],
    log: "logs/miniCon/ava2clust/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/miniCon/ava2clust/{barcode}_{c1}_{c2}_{c3}.txt"
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell:
        "python {workflow.basedir}/scripts/miniclust.py -p {params.prefix}"
        " -R {params.max_recurs}"
        " -s {params.min_score_frac} -n {params.min_reads} {input} > {log} 2>& 1"

checkpoint cls_miniCon:
    input:
        "clust/clusters", ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: get_kmerclust(wc),
        clss = lambda wc: expand("miniCon/ava2clust/{bc_cls}.csv", bc_cls = glob_wildcards(checkpoints.cls_meshclust.get(**wc).output[0] + "/{bc_cls}.csv").bc_cls),
    output: directory("miniCon/clusters")
    params:
        min_size = config["min_support_reads"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.clss):
            num_lines = sum(1 for line in open(i))
            if num_lines > params.min_size:
                bc_clss = i.split("/")[-1].removesuffix(".csv")
                df = pd.read_csv(i)
                # cluster starts from 0
                #df["cluster"] = df["cluster"] - 1
                for cand, df_cand in df.groupby('cluster'):
                    if len(df_cand) >= params.min_size:
                        bc_clss_cand = "/{bc_clss}cand{cand}".format(bc_clss=bc_clss, cand=cand)
                        df_cand['read_id'].to_csv(output[0] + bc_clss_cand + ".csv", header = False, index = False)
                        centroid_idx  = df_cand['clust_read_score'].idxmax()
                        centroid_read = df_cand.loc[centroid_idx, ['read_id']]
                        centroid_read.to_csv(output[0] + bc_clss_cand + ".centroid", header=False, index=False)
   
use rule fqs_split_meshclust as fqs_split_miniCon with:
    input:
        members = "miniCon/clusters/{barcode}_{c1}_{c2}_{c3}cand{cand}.csv",
        centroid = "miniCon/clusters/{barcode}_{c1}_{c2}_{c3}cand{cand}.centroid",
        split = rules.fqs_split_meshclust.output.members,
    output:
        members = temp("miniCon/split/{barcode}_{c1}_{c2}_{c3}cand{cand}.fastq"),
        centroid = temp("miniCon/polish/{barcode}_{c1}_{c2}_{c3}cand{cand}/minimap2/raw.fna"),

# isoCon
def check_isoCon_batch(batch_size = config["IsoCon"]["max_batch_size"]):
    check_val("IsoCon batch_size", batch_size, int)
    if int(batch_size) < -1:
        raise ValueError("IsoCon batch_size only accepts integer >= -1.")
check_isoCon_batch()

rule run_isoCon:
    input: rules.fqs_split_meshclust.output.members
    output:
        cls = temp("isoCon/{barcode}_{c1}_{c2}_{c3}/IsoCon/cluster_info.tsv"),
        fna = temp("isoCon/{barcode}_{c1}_{c2}_{c3}/IsoCon/final_candidates.fa"),
    conda: "../envs/isONcorCon.yaml"
    params:
        prefix = "isoCon/{barcode}_{c1}_{c2}_{c3}",
        min_candidates = config["min_support_reads"],
        neighbor_search_depth =  int(config["IsoCon"]["neighbor_search_depth"]) if config["IsoCon"]["neighbor_search_depth"] else 2**32,
        p_value_threshold = config["IsoCon"]["p_value_threshold"],
        max_batch_size = -1 if int(config["IsoCon"]["max_batch_size"]) == -1 else int(config["IsoCon"]["max_batch_size"]) * 4,
    log: "logs/isoCon/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/isoCon/{barcode}_{c1}_{c2}_{c3}.txt"
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
          split -l $nlines_per_batch -a3 -d --additional-suffix='.fastq' {input} {params.prefix}/IsoCon/batches/b >{log} 2>&1
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

def get_isoCon(wildcards):
    bc_clss = glob_wildcards(checkpoints.cls_meshclust.get(**wildcards).output[0] + "/{bc_cls}.csv").bc_cls
    return {
    "clss": expand("isoCon/{bc_cls}/IsoCon/cluster_info.tsv", bc_cls = bc_clss),
    "cands": expand("isoCon/{bc_cls}/IsoCon/final_candidates.fa", bc_cls = bc_clss)
    }

checkpoint cls_isoCon:
    input: 
        "clust/clusters", ".qc_DONE",
        lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc)),
        lambda wc: get_kmerclust(wc),
        unpack(get_isoCon),
    output: directory("isoCon/clusters"),
    params:
        min_candidates = config["min_support_reads"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.clss):
            # if empty, skip
            if os.stat(i).st_size == 0:
                continue
            bc_clss = i.split('/')[-3]
            df_i = pd.read_csv(i, sep = '\t', header = None, usecols=range(2))
            df_i.columns = ['seqid', 'cluster']
            # only leave candidate id 'transcript_id_support_num'
            df_i['cluster'] = df_i['cluster'].apply(lambda x: x.split('_')[1])
            for cand, df_clust in df_i.groupby('cluster'):
                if len(df_clust) < params.min_candidates:
                    continue
                df_clust['seqid'].to_csv(output[0] + "/{bc_clss}cand{cand}.csv".format(bc_clss=bc_clss, cand=cand), header = False, index = False)

rule get_IsoCon_cand:
    input:
        cls = "isoCon/clusters/{barcode}_{c1}_{c2}_{c3}cand{cand}.csv",
        cands = rules.run_isoCon.output.fna,
    output: temp("isoCon/polish/{barcode}_{c1}_{c2}_{c3}cand{cand}/minimap2/raw.fna")
    run:
            outdir = os.path.dirname(output[0])
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            bc_cls, cand = input.cls.removesuffix('.csv').split('/')[-1].split('cand')
            with open (input.cands, 'r') as fi:
                lines = fi.readlines()

            for i, line in enumerate(lines):
                if line.startswith('>transcript_' + cand + '_support'):
                    header = '>' + bc_cls + 'cand' + cand + '\n'
                    with open (output[0], 'w') as fo:
                        fo.write(header)
                        fo.write(lines[i+1])
                    break

use rule fqs_split_isONclust as fqs_split_isoCon with:
    input:
        cluster = "isoCon/clusters/{barcode}_{c1}_{c2}_{c3}cand{cand}.csv",
        fqs = rules.fqs_split_meshclust.output.members,
    output:
        temp("isoCon/split/{barcode}_{c1}_{c2}_{c3}cand{cand}.fastq")

# polish with racon and medaka
# reused in racon iterations
rule minimap2polish:
    input: 
      ref = "{consensus}/polish/{bc_cls_cand}/minimap2/{assembly}.fna",
      fastq = "{consensus}/split/{bc_cls_cand}.fastq",
    output: temp("{consensus}/polish/{bc_cls_cand}/minimap2/{assembly}.paf"),
    params:
        x = config["minimap2"]["x_map"]
    conda: "../envs/minimap2.yaml"
    log: "logs/{consensus}/{bc_cls_cand}/minimap2_{assembly}.log"
    benchmark: "benchmarks/{consensus}/{bc_cls_cand}/minimap2_{assembly}.txt"
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
        prefix = "{consensus}/polish/{bc_cls_cand}/minimap2/raw"
    else:
        prefix = "{consensus}/polish/{bc_cls_cand}/minimap2/racon_{iter}".format(consensus=wildcards.consensus,
            bc_cls_cand=wildcards.bc_cls_cand, iter=str(int(wildcards.iter) - 1))
    return(prefix + ".paf", prefix + ".fna")

rule racon:
    input:
        "{consensus}/split/{bc_cls_cand}.fastq",
        lambda wc: get_racon_input(wc),
    output: temp("{consensus}/polish/{bc_cls_cand}/minimap2/racon_{iter}.fna")
    params:
        m = config["racon"]["m"],
        x = config["racon"]["x"],
        g = config["racon"]["g"],
        w = config["racon"]["w"],
    conda: "../envs/racon.yaml"
    log: "logs/{consensus}/{bc_cls_cand}/racon_{iter}.log"
    benchmark: "benchmarks/{consensus}/{bc_cls_cand}/racon_{iter}.txt"
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
            fna = "{consensus}/polish/{bc_cls_cand}/minimap2/raw.fna".format(
                consensus=wildcards.consensus, bc_cls_cand=wildcards.bc_cls_cand)
        else:
            fna = "{consensus}/polish/{bc_cls_cand}/minimap2/racon_{iter}.fna".format(
                consensus=wildcards.consensus, bc_cls_cand=wildcards.bc_cls_cand, iter=str(int(racon_iter)))
    else:
        fna = "{consensus}/polish/{bc_cls_cand}/medaka_{iter}/consensus.fasta".format(
            consensus=wildcards.consensus, bc_cls_cand=wildcards.bc_cls_cand, iter=str(int(wildcards.iter2) - 1))
    if index == True:
        return(fna + ".fai", fna + ".map-ont.mmi")
    else:
        return(fna)

rule medaka_consensus:
    input:
        fna = lambda wc: get_medaka_files(wc),
        fastq = "{consensus}/split/{bc_cls_cand}.fastq",
    output: 
        temp(expand("{{consensus}}/polish/{{bc_cls_cand}}/medaka_{{iter2}}/consensus{ext}",
        ext = [".fasta", ".fasta.gaps_in_draft_coords.bed", "_probs.hdf"])),
        temp(expand("{{consensus}}/polish/{{bc_cls_cand}}/medaka_{{iter2}}/calls{ext}",
        ext = ["_to_draft.bam", "_to_draft.bam.bai"])),
    params:
        m = config["medaka"]["m"],
        cudnn = 'CUDA_VISIBLE_DEVICES=""' if config["medaka"]["cudnn"] is False else 'TF_FORCE_GPU_ALLOW_GROWTH=true',
        _dir = "{consensus}/polish/{bc_cls_cand}/medaka_{iter2}",
        inedxs = lambda wc: get_medaka_files(wc, index = True),
    conda: "../envs/medaka.yaml"
    log: "logs/{consensus}/{bc_cls_cand}/medaka_{iter2}.log"
    benchmark: "benchmarks/{consensus}/{bc_cls_cand}/medaka_{iter2}.txt"
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
            export OLD_LD_LIBRARY_PATH=${{LD_LIBRARY_PATH}}
            export LD_LIBRARY_PATH="$CONDA_PREFIX/lib":${{LD_LIBRARY_PATH}}
            export TF_CPP_MIN_LOG_LEVEL='2'

            export {params.cudnn}
            medaka_consensus -i {input.fastq} -d {input.fna} -o {params._dir} -t {threads} -m {params.m} > {log} 2>&1
            rm -f {params.inedxs}

            export LD_LIBRARY_PATH=${{OLD_LD_LIBRARY_PATH}}
            unset OLD_LD_LIBRARY_PATH
        fi
        """

# get polished asembly
def merge_consensus(fi, fo):
    with open(fo, "w") as out:
        for i in fi:
            # if fi is empty, skip; rm dummy output
            if os.stat(i).st_size == 0:
                continue
            bc_cls_cand = i.split("/")[-3]
            if "cand" not in bc_cls_cand:
                bc_cls_cand = bc_cls_cand + "cand1"
            
            with open(i, "r") as inp:
                for line in inp:
                    if line.startswith(">"):
                        line = ">" + bc_cls_cand + "\n"
                    out.write(line)

def get_candidates(wildcards):
    if wildcards.consensus == "kmerCon":
        candidates = glob_wildcards(checkpoints.cls_kmerCon.get(**wildcards).output[0] + "/{bc_cls_cand}.csv").bc_cls_cand
    elif wildcards.consensus == "miniCon":
        candidates = glob_wildcards(checkpoints.cls_miniCon.get(**wildcards).output[0] + "/{bc_cls_cand}.csv").bc_cls_cand
    elif wildcards.consensus == "isoCon":
        candidates = glob_wildcards(checkpoints.cls_isoCon.get(**wildcards).output[0] + "/{bc_cls_cand}.csv").bc_cls_cand
    elif wildcards.consensus == "umiCon":
        candidates = glob_wildcards(checkpoints.cls_umiCon.get(**wildcards).output[0] + "/{bc_cls_cand}.csv").bc_cls_cand
    else:
        raise ValueError("Unknown consensus method: " + wildcards.consensus)
    return candidates
 
def get_consensus(wildcards, medaka_iter = config["medaka"]["iter"], racon_iter = config["racon"]["iter"]):
    # iter >= 0, integer
    check_val("racon iter", racon_iter, int)
    check_val("medaka iter", medaka_iter, int)
    if racon_iter < 0 or medaka_iter < 0:
        raise ValueError("racon and medaka iter shall be >= 0")
    candidates = get_candidates(wildcards)

    if medaka_iter == 0:
        if racon_iter == 0:
            return expand("{{consensus}}/polish/{cand}/minimap2/raw.fna", consensus = wildcards.consensus, cand = candidates)
        else:
            return expand("{{consensus}}/polish/{cand}/minimap2/racon_{iter}.fna", consensus = wildcards.consensus, cand = candidates, iter = racon_iter)
    else:
        return expand("{{consensus}}/polish/{cand}/medaka_{iter}/consensus.fasta", consensus = wildcards.consensus, cand = candidates, iter = medaka_iter)
            
rule collect_consensus:
    input: 
        "{consensus}/clusters",
        fna = lambda wc: get_consensus(wc),
    output: "{consensus}/{consensus}.fna"
    run: merge_consensus(fi = input.fna, fo = output[0])
