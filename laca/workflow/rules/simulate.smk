# expand list from comma separated string
def exp_list(x):
    return [y.strip() for y in x.split(",")]
minids = exp_list(config["simulate"]["min_id"])
maxids = exp_list(config["simulate"]["max_id"])
quantities = exp_list(config["simulate"]["badread"]["quantity"])

depths = exp_list(config["simulate"]["depth"])
def check_depth(n=config["simulate"]["n_refs"], depth=depths):
    if len(depth) != n:
        raise ValueError("Number of depth values must be equal to number of references")
check_depth()

# link to download database
dict_db = {
    #"greengene": "https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_tax_with_hOTUs_99_reps.fasta",
    "silva": "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
}

def get_val(key, dict_i):
    return dict_i[key]

rule download_markerDB:
    output: temp("simulate/{db}/rep.fna")
    params:
        _dir = lambda wc: os.path.join(os.getcwd() + "/simulate", wc.db),
        db_links = lambda wc: get_val(wc.db, dict_db),
    log: "logs/simulate/{db}/download_markerDB.log"
    benchmark: "benchmarks/simulate/{db}/download_markerDB.txt"
    shell:
        """
        mkdir -p {params._dir}
        wget -P {params._dir} {params.db_links} > {log} 2>&1
        
        file=$(ls {params._dir})
        if [[ $file == *.gz ]]
        then
            gunzip {params._dir}/$file
            file=$(ls {params._dir})
        fi
        mv {params._dir}/$file {output}
        """
# exclude sequences with "N", rna2dna
rule rna2dna:
    input: rules.download_markerDB.output
    output: "simulate/{db}/repDNA.fna"
    conda: "../envs/seqkit.yaml"
    log: "logs/simulate/{db}/rna2dna.log"
    benchmark: "benchmarks/simulate/{db}/rna2dna.txt"
    shell: "seqkit grep -s -v -i -p N {input} | seqkit seq --rna2dna -w0 -o {output} 2> {log}"

# two round clustering to select refs
# first round: pick seqs in the same cluster, with min identity
# second round: pick representatives of each clusters, with max identity
rule cls_ref1:
    input: rules.rna2dna.output
    output:
        rep = temp("simulate/{db}/cls_ref/id_{minid}/mmseqs2_rep_seq.fasta"),
        all_by_cluster = temp("simulate/{db}/cls_ref/id_{minid}/mmseqs2_all_seqs.fasta"),
        tsv = temp("simulate/{db}/cls_ref/id_{minid}/mmseqs2_cluster.tsv"),
        tmp = temp(directory("simulate/tmp/{db}/cls_ref/id_{minid}")),
    conda: "../envs/mmseqs2.yaml"
    params:
        minid = lambda wc: int(wc.minid) * 0.01,
        prefix = lambda wc: "simulate/{db}/cls_ref/id_{minid}/mmseqs2".format(db=wc.db, minid=wc.minid),
        c = 1,
    log: "logs/simulate/{db}/cls_ref/id_{minid}.log"
    benchmark: "benchmarks/simulate/{db}/cls_ref/id_{minid}.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input} "
        "{params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.minid} -c {params.c} --cluster-reassign > {log} 2>&1"

# take the largest cluster for seconda clustering (min id)
rule cls_pick:
    input: rules.cls_ref1.output.all_by_cluster
    output: temp("simulate/{db}/cls_ref/id_{minid}/mmseqs2_cluster_top.fasta")
    run:
       # fasta: two lines per sequence: header with ">" and sequence;
       # only one line header suggests a new cluster append to dictory
       with open(input[0], "r") as f:
            d = {}
            header = []
            for line in f:
                if line.startswith(">"):
                    line1 = line
                    line2 = next(f)
                    if line2.startswith(">"):
                        header.append(line1)
                        d[header[-1]] = line2 + next(f)
                    else:
                        d[header[-1]] += line1 + line2
                 
            # get the largest cluster， "/n" in the values and return the key
            max_cluster = max(d.values(), key=lambda x: x.count("\n"))           
            # get the header of the largest cluster
            max_cluster_header = [k for k, v in d.items() if v == max_cluster]
            # write to output
            with open(output[0], "w") as f2:
                f2.write(d[max_cluster_header[0]])
                
# take representatives from 2nd clustering, (max id)
use rule cls_ref1 as cls_ref2 with:
    input: rules.cls_pick.output
    output: 
        rep = "simulate/{db}/cls_ref/id_{minid}_{maxid}/mmseqs2_rep_seq.fasta",
        all_by_cluster = temp("simulate/{db}/cls_ref/id_{minid}_{maxid}/mmseqs2_all_seqs.fasta"),
        tsv = temp("simulate/{db}/cls_ref/id_{minid}_{maxid}/mmseqs2_cluster.tsv"),
        tmp = temp(directory("simulate/tmp/{db}/cls_ref/id_{minid}_{maxid}")),
    params:
        minid = lambda wc: int(wc.maxid) * 0.01,
        prefix = lambda wc: "simulate/{db}/cls_ref/id_{minid}_{maxid}/mmseqs2".format(db = wc.db, minid=wc.minid, maxid=wc.maxid),
        c = 1,
    log: "logs/simulate/{db}/cls_ref/id_{minid}_{maxid}.log"
    benchmark: "benchmarks/simulate/{db}/cls_ref/id_{minid}_{maxid}.txt"

# set a length range
rule filter_ref_len:
    input: rules.cls_ref2.output.rep    
    output: temp("simulate/{db}/cls_ref/id_{minid}_{maxid}/mmseqs2_rep_seq_f.fasta")
    conda: "../envs/seqkit.yaml"
    params:
        m = config["simulate"]["min_len"],
        M = config["simulate"]["max_len"],
    log: "logs/simulate/{db}/cls_ref/id_{minid}_{maxid}_filter_len.log"  
    benchmark: "benchmarks/simulate/{db}/cls_ref/id_{minid}_{maxid}_filter_len.txt"
    threads: config["threads"]["normal"]
    shell: "cat {input} | seqkit seq -j {threads} -m {params.m} -M {params.M} -i -w0 > {output} 2> {log}"

# subsample ref
rule subsample_cls_ref:
    input: rules.filter_ref_len.output
    output:
        p = temp("simulate/{db}/cls_ref/id_{minid}_{maxid}/subsampled_p.fasta"),
        n = temp("simulate/{db}/cls_ref/id_{minid}_{maxid}/subsampled.fasta"),
    conda: "../envs/seqkit.yaml"
    params:
        p = 1,
        n = config["simulate"]["n_refs"],
    log: "logs/simulate/{db}/cls_ref/id_{minid}_{maxid}_subsample.log"
    benchmark: "benchmarks/simulate/{db}/cls_ref/id_{minid}_{maxid}_subsample.txt"
    threads: 1
    shell:
        """
        seqkit sample -p {params.p} -j {threads} {input} -w0 -o {output.p} -s123 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -w0 -o {output.n} 2>> {log}
        """

# keep fasta header unique
rule reheader:
    input: rules.subsample_cls_ref.output.n
    output: "simulate/{db}/cls_ref/id_{minid}_{maxid}/ref.fasta"
    params:
        depth = depths,
    log: "logs/simulate/{db}/cls_ref/id_{minid}_{maxid}_reheader.log"
    benchmark: "benchmarks/simulate/{db}/cls_ref/id_{minid}_{maxid}_reheader.txt"
    run:
        with open(output[0], "w") as out:
            with open (input[0], "r") as inp:
                i = 1
                for line in inp:
                    if line.startswith(">"):
                        line = ">ref_{i} depth={depth_i}\n".format(i=i, depth_i=params.depth[i-1])
                        i += 1 
                    out.write(line)
        # if i < len(params.depth), only depth[:i] is appended and raise warning to log
        if i < len(params.depth):
            with open(log[0], "w") as logfile:
                logfile.write(
                    "Number of references is less than number of depth values,\nonly the first {i} depth values in config are appended to fasta header.".format(i=i))

def calc_mean_and_stdev(fasta):
    lengths = []
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            lengths.append(len(line.strip()))

    mean = sum(lengths) / len(lengths)
    variance = sum((x - mean) ** 2 for x in lengths) / len(lengths)
    stdev = variance ** 0.5
    
    return "{mean},{stdev}".format(mean=mean, stdev=stdev)

# simulate long reads with badread
rule badread:
    input: rules.reheader.output
    output: "simulate/{db}/badread/id_{minid}_{maxid}/reads{quantity}.fastq"
    conda: "../envs/badread.yaml"
    params:
        length = lambda wc, input: calc_mean_and_stdev(fasta=input[0]),
        identity = config["simulate"]["badread"]["identity"],
        error_model = config["simulate"]["badread"]["error_model"],
        qscore_model = config["simulate"]["badread"]["qscore_model"],
        seed = 123,
        start_adapter = config["simulate"]["badread"]["start_adapter"],
        end_adapter = config["simulate"]["badread"]["end_adapter"],
        start_adapter_seq = config["simulate"]["badread"]["start_adapter_seq"],
        end_adapter_seq = config["simulate"]["badread"]["end_adapter_seq"],
        junk_reads = config["simulate"]["badread"]["junk_reads"],
        random_reads = config["simulate"]["badread"]["random_reads"],
        chimeras = config["simulate"]["badread"]["chimeras"],
        glitches = config["simulate"]["badread"]["glitches"],
    log: "logs/simulate/{db}/badread/id_{minid}_{maxid}/reads{quantity}.log"
    benchmark: "benchmarks/simulate/{db}/badread/id_{minid}_{maxid}/reads{quantity}.txt"
    shell:
        """
        badread simulate --reference {input} --quantity {wildcards.quantity} \
            --length {params.length} --identity {params.identity} \
            --error_model {params.error_model} --qscore_model {params.qscore_model} \
            --seed {params.seed} --start_adapter {params.start_adapter} \
            --end_adapter {params.end_adapter} --start_adapter_seq {params.start_adapter_seq} \
            --end_adapter_seq {params.end_adapter_seq} --junk_reads {params.junk_reads} \
            --random_reads {params.random_reads} --chimeras {params.chimeras} \
            --glitches {params.glitches} > {output} 2> {log}
        """

localrules: simulate
rule simulate:
    input: expand("simulate/{db}/badread/id_{mixid}/reads{quantity}.fastq", db = [k for k in dict_db.keys()], mixid = [f"{minid}_{maxid}" for (minid, maxid) in zip(minids, maxids)], quantity=quantities)
    output: directory("demultiplexed"),
    run:
        # if demultiplexed directory exists, exit with error
        if os.path.exists(output[0]):
            raise Exception("  Demultiplexed directory ({demux_dir}) already exists.\n  To avoid unwanted overwriting, please manually remove it and re-run 'laca simulate'.".format(demux_dir=output[0]))
        else:
            os.makedirs(output[0])
        for fq in input:
            # no wildcards, use input path, keep digits only
            # last second directory name, trim "id_", replace "_" with ""
            mixid_digits = os.path.basename(os.path.dirname(fq)).replace("id_", "").replace("_", "")
            # file name, trim ".fastq", "reads"
            quantity = os.path.basename(fq).removesuffix(".fastq").replace("reads", "")
            quantity_digits = ''.join(str(i) for i in quantity if i.isdigit())
            dirname = os.path.join(output[0], "si" + mixid_digits + quantity_digits) 
            os.makedirs(dirname)
            with open (fq, "r") as fi:
                with open (os.path.join(dirname, "badread.fastq"), "w") as fo:
                    for line in fi:
                        fo.write(line)        