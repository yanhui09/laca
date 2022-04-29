# Prepare nanopore reads
REF = config["nanosim"]["ref"]
# expand list from comma separated string
def exp_list(x):
    return [y.strip() for y in x.split(",")]
IDS = exp_list(config["nanosim"]["id"])
NS = exp_list(config["nanosim"]["simulator"]["n"])

# characterize error model
NREAD = config["nanosim"]["read_analysis"]["read"]
RG = config["nanosim"]["read_analysis"]["rg"]

rule read_analysis:
    input: 
      read = NREAD,
      ref = RG,
    output: OUTPUT_DIR + "/nanosim/model/nanosim_model_profile"
    conda: "../envs/nanosim.yaml"
    params:
        a = config["nanosim"]["read_analysis"]["a"],
        prefix = OUTPUT_DIR + "/nanosim/model/nanosim",
    log: OUTPUT_DIR + "/logs/nanosim/read_analysis.log"
    benchmark: OUTPUT_DIR + "/benchmarks/nanosim/read_analysis.txt"
    threads: config["nanosim"]["threads"]
    shell: "read_analysis.py genome -i {input.read} -rg {input.ref} -o {output} -a {params.a} -t {threads} > {log} 2>&1"

# cluster ref by identity
rule cls_ref:
    input: REF
    output: 
        rep = OUTPUT_DIR + "/nanosim/cls_ref/id_{minid}/mmseqs_rep_seq.fasta",
        all_by_cluster = temp(OUTPUT_DIR + "/nanosim/cls_ref/id_{minid}/mmseqs_all_seqs.fasta"),
        tsv = temp(OUTPUT_DIR + "/nanosim/cls_ref/id_{minid}/mmseqs_cluster.tsv"),
        tmp = temp(directory(OUTPUT_DIR + "/tmp/cls_ref/id_{minid}")),
    conda: "../envs/mmseqs2.yaml"
    params:
        minid = lambda wc: int(wc.minid) * 0.01,
        prefix = lambda wc: OUTPUT_DIR + "/nanosim/cls_ref/id_{minid}/mmseqs".format(minid=wc.minid),
        c = 1,
    log: OUTPUT_DIR + "/logs/nanosim/id_{minid}/cls_ref.log"
    benchmark: OUTPUT_DIR + "/benchmarks/nanosim/id_{minid}/cls_ref.txt"
    threads: config["nanosim"]["threads"]
    shell:
        "mmseqs easy-cluster {input} {params.prefix} {output.tmp} "
        "--threads {threads} --min-seq-id {params.minid} -c {params.c} > {log} 2>&1"

# subsample clustered ref
rule subsample_cls_ref:
    input: rules.cls_ref.output.rep
    output:
        p = temp(OUTPUT_DIR + "/nanosim/cls_ref/id_{minid}/subsampled_p.fasta"),
        n = temp(OUTPUT_DIR + "/nanosim/cls_ref/id_{minid}/subsampled.fasta"),
    conda: "../envs/seqkit.yaml"
    params:
        p = 1,
        n = config["nanosim"]["subsample_n"],
    log: OUTPUT_DIR + "/logs/nanosim/id_{minid}/subsample_cls_ref.log"
    benchmark: OUTPUT_DIR + "/benchmarks/nanosim/id_{minid}/subsample_cls_ref.txt"
    threads: 1
    shell:
        """
        seqkit sample -p {params.p} -j {threads} {input} -o {output.p} 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -o {output.n} 2>> {log}
        """

# keep fasta header unique
rule reheader:
    input: rules.subsample_cls_ref.output
    output: OUTPUT_DIR + "/nanosim/cls_ref/id_{minid}/ref.fasta"
    log: OUTPUT_DIR + "/logs/nanosim/cls_ref/id_{minid}/reheader.log"
    benchmark: OUTPUT_DIR + "/benchmarks/nanosim/cls_ref/id_{minid}/reheader.txt"
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
        ref = rules.reheader.output,
    output: 
        expand(
            OUTPUT_DIR + "/nanosim/simulate/{{minid}}_{{n}}/simulated{suffix}",
            suffix=["_aligned_error_profile", "_aligned_reads.fasta", "_unaligned_reads.fasta"]
        )
    conda: "../envs/nanosim.yaml"
    params:
        c = OUTPUT_DIR + "/nanosim/model/nanosim",
        o = lambda wc: OUTPUT_DIR + "/nanosim/simulate/{minid}_{n}/simulated".format(minid=wc.minid, n=wc.n),
        n = lambda wc: int(wc.n),
    log: OUTPUT_DIR + "/logs/nanosim/read_simulate/{minid}_{n}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/nanosim/read_simulate/{minid}_{n}.txt"
    threads: config["threads"]["large"]
    shell:
        "simulator.py genome -rg {input.ref} -c {params.c} -o {params.o}"
        " -n {params.n} -t {threads} -dna_type linear"
        " > {log} 2>&1"

rule nanosim:
    input: expand(OUTPUT_DIR + "/nanosim/simulate/{minid}_{n}/simulated_aligned_reads.fasta", minid=IDS, n=NS)
    output: temp(touch(OUTPUT_DIR + ".simulated_DONE"))
