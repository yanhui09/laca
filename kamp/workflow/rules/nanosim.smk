# Nanosim seems not to properly simulate the ends of reference
# https://github.com/bcgsc/NanoSim/issues/101

# expand list from comma separated string
def exp_list(x):
    return [y.strip() for y in x.split(",")]
IDS = exp_list(config["nanosim"]["id"])
NS = exp_list(config["nanosim"]["simulator"]["n"])

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

# cluster ref by identity
rule cls_ref:
    output: 
        rep = temp("nanosim/cls_ref/id_{minid}/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp("nanosim/cls_ref/id_{minid}/mmseqs_all_seqs.fasta"),
        tsv = temp("nanosim/cls_ref/id_{minid}/mmseqs_cluster.tsv"),
        tmp = temp(directory("nanosim/tmp/cls_ref/id_{minid}")),
    conda: "../envs/mmseqs2.yaml"
    params:
        minid = lambda wc: int(wc.minid) * 0.01,
        prefix = lambda wc: "nanosim/cls_ref/id_{minid}/mmseqs".format(minid=wc.minid),
        c = 1,
    log: "logs/nanosim/id_{minid}/cls_ref.log"
    benchmark: "benchmarks/nanosim/id_{minid}/cls_ref.txt"
    threads: config["threads"]["large"]
    shell:
        "mmseqs easy-cluster {workflow.basedir}/../resources/data/silva_id100subs1000.fna "
        "{params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.minid} -c {params.c} > {log} 2>&1"

# subsample clustered ref
rule subsample_cls_ref:
    input: rules.cls_ref.output.rep
    output:
        p = temp("nanosim/cls_ref/id_{minid}/subsampled_p.fasta"),
        n = temp("nanosim/cls_ref/id_{minid}/subsampled.fasta"),
    conda: "../envs/seqkit.yaml"
    params:
        p = 1,
        n = config["nanosim"]["subsample_n"],
    log: "logs/nanosim/id_{minid}/subsample_cls_ref.log"
    benchmark: "benchmarks/nanosim/id_{minid}/subsample_cls_ref.txt"
    threads: 1
    shell:
        """
        seqkit sample -p {params.p} -j {threads} {input} -o {output.p} 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -o {output.n} 2>> {log}
        """

# keep fasta header unique
rule reheader:
    input: rules.subsample_cls_ref.output.n
    output: "nanosim/cls_ref/id_{minid}/ref.fasta"
    log: "logs/nanosim/cls_ref/id_{minid}/reheader.log"
    benchmark: "benchmarks/nanosim/cls_ref/id_{minid}/reheader.txt"
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
            "nanosim/simulate/{{minid}}_{{n}}/simulated{suffix}",
            suffix=["_aligned_error_profile", "_aligned_reads.fastq", "_unaligned_reads.fastq"]
        )
    conda: "../envs/nanosim.yaml"
    params:
        o = lambda wc: "nanosim/simulate/{minid}_{n}/simulated".format(minid=wc.minid, n=wc.n),
        n = lambda wc: int(wc.n),
    log: "logs/nanosim/read_simulate/{minid}_{n}.log"
    benchmark: "benchmarks/nanosim/read_simulate/{minid}_{n}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        MODEL={input[0]}
        c=${{MODEL//_model_profile/}}
        simulator.py genome -rg {input.ref} -c "$c" -o {params.o} -n {params.n} -b guppy --fastq -t {threads} -dna_type linear --seed 1234 > {log} 2>&1
        """

def sim_demult_flag(demult="guppy"):
    # if demult != "guppy" | "minibar" | "nanosim", raise value error
    if demult != "guppy" and demult != "minibar":
        raise ValueError("Demultiplexer not recognized. Choose guppy or minibar in config.")
    return "." + demult + "_DONE"

# pseudo demultiplex
rule nanosim:
    input: expand("nanosim/simulate/{minid}_{n}/simulated_aligned_reads.fastq", minid=IDS, n=NS)
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

