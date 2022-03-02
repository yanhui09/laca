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
    output: directory(OUTPUT_DIR + "/nanosim/")
    conda: "../envs/nanosim.yaml"
    params:
        a = config["nanosim"]["read_analysis"]["a"],
    log: OUTPUT_DIR + "/logs/nanosim/read_analysis.log"
    benchmark: OUTPUT_DIR + "/benchmarks/nanosim/read_analysis.txt"
    threads: config["nanosim"]["threads"]
    shell: "read_analysis.py genome -i {input.read} -rg {input.ref} -o {output} -a {params.a} -t {threads} > {log} 2>&1"

