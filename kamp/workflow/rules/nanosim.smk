# Nanosim seems not to properly simulate the ends of reference
# https://github.com/bcgsc/NanoSim/issues/101

# expand list from comma separated string
def exp_list(x):
    return [y.strip() for y in x.split(",")]
minids = exp_list(config["nanosim"]["min_seq_id"])
maxids = exp_list(config["nanosim"]["max_seq_id"])
ns = exp_list(config["nanosim"]["simulator"]["n"])

rule read_analysis:
    output: "nanosim/model/training_model_profile"
    conda: "../envs/nanosim.yaml"
    params:
        prefix = os.path.join(os.getcwd(), "nanosim/model/training"),
    log: "logs/nanosim/read_analysis.log"
    benchmark: "benchmarks/nanosim/read_analysis.txt"
    threads: config["threads"]["large"]
    shell: 
        "read_analysis.py genome -i {workflow.basedir}/../resources/data/raw.fastq.gz "
        "-rg {workflow.basedir}/../resources/data/zymock_16s.fna "
        "-o {params.prefix} -a minimap2 -t {threads} "
        "> {log} 2>&1 "

# link to download database
dict_db = {
    #"greengene": "https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_tax_with_hOTUs_99_reps.fasta",
    "silva": "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
}

def get_val(key, dict_i):
    return dict_i[key]

rule download_markerDB:
    output: temp("nanosim/{db}/rep.fna")
    params:
        _dir = lambda wc: os.path.join(os.getcwd() + "/nanosim", wc.db),
        db_links = lambda wc: get_val(wc.db, dict_db),
    log: "logs/nanosim/{db}/download_markerDB.log"
    benchmark: "benchmarks/nanosim/{db}/download_markerDB.txt"
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

rule rna2dna:
    input: rules.download_markerDB.output
    output: "nanosim/{db}/repDNA.fna"
    conda: "../envs/seqkit.yaml"
    log: "logs/nanosim/{db}/rna2dna.log"
    benchmark: "benchmarks/nanosim/{db}/rna2dna.txt"
    shell: "seqkit seq -w0 --rna2dna {input} > {output} 2> {log}"

# two round clustering to select refs
# first round: pick seqs in the same cluster, with min identity
# second round: pick representatives of each clusters, with max identity
rule cls_ref1:
    input: rules.rna2dna.output
    output:
        rep = temp("nanosim/{db}/cls_ref/id_{minid}/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp("nanosim/{db}/cls_ref/id_{minid}/mmseqs_all_seqs.fasta"),
        tsv = temp("nanosim/{db}/cls_ref/id_{minid}/mmseqs_cluster.tsv"),
        tmp = temp(directory("nanosim/tmp/{db}/cls_ref/id_{minid}")),
    conda: "../envs/mmseqs2.yaml"
    params:
        minid = lambda wc: int(wc.minid) * 0.01,
        prefix = lambda wc: "nanosim/{db}/cls_ref/id_{minid}/mmseqs".format(db=wc.db, minid=wc.minid),
        c = 1,
    log: "logs/nanosim/{db}/cls_ref/id_{minid}.log"
    benchmark: "benchmarks/nanosim/{db}/cls_ref/id_{minid}.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {input} "
        "{params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.minid} -c {params.c} --cluster-reassign > {log} 2>&1"

# take the largest cluster for seconda clustering (min id)
rule cls_pick:
    input: rules.cls_ref1.output.all_by_cluster
    output: temp("nanosim/{db}/cls_ref/id_{minid}/mmseqs_cluster_top.fasta")
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
        rep = temp("nanosim/{db}/cls_ref/id_{minid}_{maxid}/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp("nanosim/{db}/cls_ref/id_{minid}_{maxid}/mmseqs_all_seqs.fasta"),
        tsv = temp("nanosim/{db}/cls_ref/id_{minid}_{maxid}/mmseqs_cluster.tsv"),
        tmp = temp(directory("nanosim/tmp/{db}/cls_ref/id_{minid}_{maxid}")),
    params:
        minid = lambda wc: int(wc.maxid) * 0.01,
        prefix = lambda wc: "nanosim/{db}/cls_ref/id_{minid}_{maxid}/mmseqs".format(db = wc.db, minid=wc.minid, maxid=wc.maxid),
        c = 1,
    log: "logs/nanosim/{db}/cls_ref/id_{minid}_{maxid}.log"
    benchmark: "benchmarks/nanosim/{db}/cls_ref/id_{minid}_{maxid}.txt"

# set a length range
rule filter_ref_len:
    input: rules.cls_ref2.output.rep    
    output: temp("nanosim/{db}/cls_ref/id_{minid}_{maxid}/mmseqs_rep_seq_f.fasta")
    conda: "../envs/seqkit.yaml"
    params:
        m = config["nanosim"]["min_len"],
        M = config["nanosim"]["max_len"],
    log: "logs/nanosim/{db}/cls_ref/id_{minid}_{maxid}_filter_len.log"  
    benchmark: "benchmarks/nanosim/{db}/cls_ref/id_{minid}_{maxid}_filter_len.txt"
    threads: config["threads"]["normal"]
    shell: "cat {input} | seqkit seq -j {threads} -m {params.m} -M {params.M} -i -w0 > {output} 2> {log}"

# subsample ref
rule subsample_cls_ref:
    input: rules.filter_ref_len.output
    output:
        p = temp("nanosim/{db}/cls_ref/id_{minid}_{maxid}/subsampled_p.fasta"),
        n = temp("nanosim/{db}/cls_ref/id_{minid}_{maxid}/subsampled.fasta"),
    conda: "../envs/seqkit.yaml"
    params:
        p = 1,
        n = config["nanosim"]["subsample_n"],
    log: "logs/nanosim/{db}/cls_ref/id_{minid}_{maxid}_subsample.log"
    benchmark: "benchmarks/nanosim/{db}/cls_ref/id_{minid}_{maxid}_subsample.txt"
    threads: 1
    shell:
        """
        seqkit sample -p {params.p} -j {threads} {input} -w0 -o {output.p} 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -w0 -o {output.n} 2>> {log}
        """

# keep fasta header unique
rule reheader:
    input: rules.subsample_cls_ref.output.n
    output: "nanosim/{db}/cls_ref/id_{minid}_{maxid}/ref.fasta"
    log: "logs/nanosim/{db}/cls_ref/id_{minid}_{maxid}_reheader.log"
    benchmark: "benchmarks/nanosim/{db}/cls_ref/id_{minid}_{maxid}_reheader.txt"
    run:
        with open(output[0], "w") as out:
            with open (input[0], "r") as inp:
                i = 1
                for line in inp:
                    if line.startswith(">"):
                        line = ">ref_" + str(i) + "\n"
                        i += 1 
                    out.write(line)

# simulate reads by number of reads
rule read_simulate:
    input:
        rules.read_analysis.output,
        ref = rules.reheader.output,
    output: 
        expand(
            "nanosim/{{db}}/simulate/{{minid}}_{{maxid}}_{{n}}/simulated{suffix}",
            suffix=["_aligned_error_profile", "_aligned_reads.fastq"]
        ),
        directory("nanosim/{db}/simulate/{minid}_{maxid}_{n}"),
    conda: "../envs/nanosim.yaml"
    params:
        o = lambda wc: "nanosim/{db}/simulate/{minid}_{maxid}_{n}/simulated".format(db=wc.db, minid=wc.minid, maxid=wc.maxid, n=wc.n),
        n = lambda wc: int(wc.n),
        min_len = int(config["nanosim"]["min_len"]) - 100,
        max_len = int(config["nanosim"]["max_len"]) + 100,
    log: "logs/nanosim/{db}/read_simulate/{minid}_{maxid}_{n}.log"
    benchmark: "benchmarks/nanosim/{db}/read_simulate/{minid}_{maxid}_{n}.txt"
    threads: config["threads"]["normal"]
    shell:
        """
        MODEL={input[0]}
        c=${{MODEL//_model_profile/}}
        simulator.py genome -rg {input.ref} -c "$c" -o {params.o} -n {params.n} -b guppy --fastq -t {threads} -dna_type linear --seed 123 -max {params.max_len} -min {params.min_len} > {log} 2>&1
        """

def sim_demult_flag(demult="guppy"):
    # if demult != "guppy" | "minibar" | "nanosim", raise value error
    if demult != "guppy" and demult != "minibar":
        raise ValueError("Demultiplexer not recognized. Choose guppy or minibar in config.")
    return "." + demult + "_DONE"

# pseudo demultiplex
rule nanosim:
    input: expand("nanosim/{db}/simulate/{mixid}_{n}/simulated_aligned_reads.fastq", db = [k for k in dict_db.keys()], mixid = [f"{minid}_{maxid}" for (minid, maxid) in zip(minids, maxids)], n=ns)
    output: 
        touch(".simulated_DONE"),
        directory("demultiplexed"),
        touch(sim_demult_flag(config["demultiplex"])),
    run:
        import shutil
        # replace fqs to follow demultiplexing wildcards
        if not os.path.exists(output[1]):
            os.makedirs(output[1])
        for fq in input:
            fname = os.path.basename(fq)
            # last second directory name
            dirname = os.path.basename(os.path.dirname(fq))
            # remove "_"
            dirname = dirname.replace("_", "")
            dirname = os.path.join(output[1], "si" + dirname)
            os.makedirs(dirname)
            shutil.copy(fq, os.path.join(dirname, fname))

