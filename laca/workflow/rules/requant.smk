# check the existence of the requant directory
def requant_check():
    requant_path = os.getcwd() + '/f2requant'
    if not os.path.exists(requant_path):
        raise ValueError("\n  f2requant directory not found.\n\tMake sure {} is loaded.\n".format(requant_path))
    else:
        fs = [f for f in os.listdir(requant_path) if os.path.isfile(os.path.join(requant_path, f))]
        # file with suffix .fasta, fa, or fna
        fas = [f for f in fs if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.fna')]
        # file with suffix .fastq. .fq
        fqs = [f for f in fs if f.endswith('.fastq')]
        if len(fas) == 0 or len(fqs) == 0:
            raise ValueError("\n  The required fasta (.fasta|.fa|.fna) or fastq (.fastq|.fq) files \n  are not found in {}.\n".format(requant_path))
        else:
            return True

checkpoint requant_dir:
    output: touch(directory("f2requant"))
    run:
        requant_check()

def col_info_requant(wildcards, type):
    if type == "fq":
        f = glob_wildcards(checkpoints.requant_dir.get(**wildcards).output[0]
        + "/{f, (.*).(fastq|fq)}").f
    elif type == "fa":
        f = glob_wildcards(checkpoints.requant_dir.get(**wildcards).output[0]
        + "/{f, (.*).(fasta|fna|fa)}").f
    else:
        raise ValueError("\n  Only fa, fq as valid types.\n")
    f = sorted(set(f))
    return expand("f2requant/{f}", f=f)

rule concatenate_ref:
    input: lambda wildcards: col_info_requant(wildcards, "fa"),
    output: temp("requant/combined.fna")
    log: "logs/requant/concatenate_ref.log"
    benchmark: "benchmarks/requant/concatenate_ref.txt"
    shell:
        "cat {input} > {output} 2> {log}"

use rule derep_denoised_seqs as derep_denoised_seqs_re with:
    input: 
        first = rules.concatenate_ref.output,
    output: 
        rep = temp("requant/mmseqs2_rep_seq.fasta"),
        all_by_cluster = temp("requant/mmseqs2_all_seqs.fasta"),
        tsv = temp("requant/mmseqs2_cluster.tsv"),
        tmp = temp(directory("tmp")),
    params:
        prefix = "requant/mmseqs2",
        mid = config["mmseqs2"]["min_id"],
        c = config["mmseqs2"]["c"],
    log: 
        "logs/requant/derep_denoised_seqs.log"
    benchmark: 
        "benchmarks/requant/derep_denoised_seqs.txt"

# rm duplicates of reverse complements
use rule rmdup_revcom as rmdup_revcom_re with:
    input: 
        rules.derep_denoised_seqs_re.output.rep
    output: 
        temp("requant/mmseqs2_rep_seq_rmdup.fasta")
    log: 
        "logs/requant/rmdup_revcom.log"
    benchmark: 
        "benchmarks/requant/rmdup_revcom.txt"

# keep fasta header unique
use rule rename_fasta_header as rename_fasta_header_re with:
    input: 
        rules.rmdup_revcom_re.output
    output: 
        "rep_seqs_requant.fasta"
    log: 
        "logs/requant/rename_fasta.log"
    benchmark: 
        "benchmarks/requant/rename_fasta.txt"

# create abundance matrix with minimap2
use rule index_repseqs as index_repseqs_re with:
    input: 
        rules.rename_fasta_header_re.output
    output: 
        mmi = temp("rep_seqs_requant.mmi"),
        dict = temp("rep_seqs_requant.dict"),
    log: 
        "logs/requant/index_repseqs.log"
    benchmark: 
        "benchmarks/requant/index_repseqs.txt"

use rule minimap2repseqs as minimap2repseqs_re with:
    input:
        fq = "f2requant/{barcode}.fastq",
        mmi = rules.index_repseqs_re.output.mmi,
        dict = rules.index_repseqs_re.output.dict,
    output: 
        bam = temp("requant/mapped/{barcode}.bam"),
        sort = temp("requant/mapped/{barcode}.sorted.bam"),
        bai = temp("requant/mapped/{barcode}.sorted.bam.bai"),
        counts = temp("requant/mapped/{barcode}.count"),
    log: 
        "logs/requant/minimap2/{barcode}.log"
    benchmark: 
        "benchmarks/requant/minimap2/{barcode}.txt"
    
def get_qout_re(wildcards, type_o):
    barcodes = col_info_requant(wildcards, "fq")
    # extract basename & strip extension
    barcodes = [os.path.basename(barcode).split(".")[0] for barcode in barcodes]
    if type_o == "bam":
        output = expand("requant/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand("requant/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand("requant/mapped/{barcode}.count", barcode=barcodes)
    else:
        raise ValueError("type_o must be 'bam', 'bai', or 'count'")
    return output

# biom format header
use rule rowname_kOTU as rowname_kOTU_re with:
    input:
        bam = lambda wildcards: get_qout_re(wildcards, "bam"),
        bai = lambda wildcards: get_qout_re(wildcards, "bai"),
    output: 
        temp("requant/rowname_seqs")
    log: 
        "logs/requant/rowname_kOTU.log"
    benchmark: 
        "benchmarks/requant/rowname_kOTU.txt"
       
use rule count_matrix as count_matrix_re with:
    input:
        rowname_seqs = rules.rowname_kOTU_re.output,
        seqs_count = lambda wildcards: get_qout_re(wildcards, "count"),
    output: 
        "count_matrix_requant.tsv"