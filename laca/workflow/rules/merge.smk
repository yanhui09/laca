# check the existence of the requant directory
def check_run(laca_wd):
    if not os.path.isdir(laca_wd):
        raise ValueError("\n  Directory not found.\n  Make sure {} exits.\n".format(laca_wd))
    else:
        # check if count_matrix.tsv and repseqs.fasta exist
        items = ["count_matrix.tsv", "rep_seqs.fasta"]
        for item in items:
            if not os.path.isfile(os.path.join(laca_wd, item)):
                raise ValueError("\n  {f} not found.\n  Do 'laca run quant' first in the working directory of {d}.\n".format(
                    f = os.path.join(laca_wd, item), d = laca_wd))

def collect_laca_runs(laca_run):
    # check if laca_run is a list
    if not isinstance(laca_run, list):
        raise ValueError("\n  'merge_runs' in config must be a YAML list.\n  Specify runs by '-' under 'merge_run' in config.\n")
        # check if laca_run is empty
        if len(laca_run) == 0:
            raise ValueError("\n  'merge_runs' in config is empty.\n  Specify runs by '-' under 'merge_run' in config.\n")
    
    seqs = []
    tables = []
    for i in laca_run:
        check_run(i)
        seqs.append(os.path.join(i, "rep_seqs.fasta"))
        tables.append(os.path.join(i, "count_matrix.tsv"))
    return {
        "seqs": seqs,
        "tables": tables
    }
    
rule merge_repseqs:
    output: temp("merge/seqs_combined.fasta")
    run:
        seqs = collect_laca_runs(laca_run = config["merge_runs"])["seqs"]
        with open(output[0], "w") as fo:
            for i in seqs:
                with open(i, "r") as fi:
                    for line in fi:
                        if line.startswith(">"):
                            # append suffix with run name, OTU_1_run_name
                            line = line.rstrip() + "_" + i.split("/")[-2] + "\n"
                        fo.write(line)
        
# dereplicate sequences with MMseqs2
use rule drep_consensus as drep_seqs_merged with:
    input: rules.merge_repseqs.output
    output: 
        rep = temp("merge/mmseqs2_rep_seq.fasta"),
        all_by_cluster = temp("merge/mmseqs2_all_seqs.fasta"),
        tsv = temp("merge/mmseqs2_cluster.tsv"),
        tmp = temp(directory("merge/tmp")),
    params:
        prefix = "merge/mmseqs2",
        mid = config["mmseqs2"]["min_id"],
        c = config["mmseqs2"]["c"],
    log: 
        "logs/merge/derep_merged_seqs.log"
    benchmark: 
        "benchmarks/merge/derep_merged_seqs.txt"
    
# rename fasta header
use rule rename_drep_seqs as rename_drep_seqs_merged with:
    input: 
        rules.drep_seqs_merged.output.rep
    output: 
       "rep_seqs_merged.fasta"

rule matrix_merged:
    input: rules.drep_seqs_merged.output.tsv
    output: "count_matrix_merged.tsv"
    run:
        import pandas as pd
        # OTU <- derep_cls -> cand_cls
        derep_cls = pd.read_csv(input[0], sep="\t", header=None)
        derep_cls.columns = ["rep_cls","cls"]
        # OTU with incremental number by rep_cls, 1, 2, ...
        derep_cls["OTU"] = derep_cls["rep_cls"].factorize()[0] + 1

        tables = collect_laca_runs(laca_run = config["merge_runs"])['tables']
        for table in list(tables):
            # read count matrix
            count_matrix = pd.read_csv(table, sep="\t", index_col=0)
            # rename columns with run name
            count_matrix.columns = [col + "_" + table.split("/")[-2] for col in count_matrix.columns]
            # 'cls' with run name appended
            count_matrix["cls"] = count_matrix.index + "_" + table.split("/")[-2]
            # merge count matrix
            derep_cls = pd.merge(derep_cls, count_matrix, on="cls", how="left")
        
        # fill NaN with 0
        derep_cls = derep_cls.fillna(0)
        # drop 'cls' and 'rep_cls', group by OTU and sum
        cls = derep_cls.drop(["cls", "rep_cls"], axis=1).groupby("OTU").sum()
        cls = cls.astype(int)
        cls.index = ["OTU_" + str(i) for i in cls.index]
        cls.index.name = "#OTU ID"
        cls.to_csv(output[0], sep="\t", header=True, index=True)
