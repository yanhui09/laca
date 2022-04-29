# check the existence of the requant directory
def requant_check(path, chimera_check = True):
    requant_path = path + '/f2requant'
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
            return chimeraF(chimera_check)

checkpoint requant_dir:
    input: requant_check(OUTPUT_DIR, config["chimeraF"])
    output: touch(directory(OUTPUT_DIR + "/f2requant"))

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
    return expand(OUTPUT_DIR + "/f2requant/{f}", f=f)

rule concatenate_ref:
    input: lambda wildcards: col_info_requant(wildcards, "fa"),
    output: temp(OUTPUT_DIR + "/requant/combined.fna")
    log: OUTPUT_DIR + "/logs/requant/concatenate_ref.log"
    benchmark: OUTPUT_DIR + "/benchmarks/requant/concatenate_ref.txt"
    shell:
        "cat {input} > {output} 2> {log}"

# dereplicate sequences with mmseqs
use rule derep_denoised_seqs as derep_denoised_seqs_re with:
    input: 
        first = rules.concatenate_ref.output,
    output: 
        rep = temp(OUTPUT_DIR + "/requant/mmseqs_rep_seq.fasta"),
        all_by_cluster = temp(OUTPUT_DIR + "/requant/mmseqs_all_seqs.fasta"),
        tsv = temp(OUTPUT_DIR + "/requant/mmseqs_cluster.tsv"),
        tmp = temp(directory(OUTPUT_DIR + "/tmp")),
    params:
        prefix = OUTPUT_DIR + "/requant/mmseqs",
        mid = config["mmseqs"]["min-seq-id"],
        c = config["mmseqs"]["c"],
    log: 
        OUTPUT_DIR + "/logs/requant/derep_denoised_seqs.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/derep_denoised_seqs.txt"

# rm duplicates of reverse complements
use rule rmdup_revcom as rmdup_revcom_re with:
    input: 
        rules.derep_denoised_seqs_re.output.rep
    output: 
        temp(OUTPUT_DIR + "/requant/mmseqs_rep_seq_rmdup.fasta")
    log: 
        OUTPUT_DIR + "/logs/requant/rmdup_revcom.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/rmdup_revcom.txt"

# keep fasta header unique
use rule rename_fasta_header as rename_fasta_header_re with:
    input: 
        rules.rmdup_revcom_re.output
    output: 
        OUTPUT_DIR + "/rep_seqs_requant.fasta"
    log: 
        OUTPUT_DIR + "/logs/requant/rename_fasta.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/rename_fasta.txt"

# create abundance matrix with minimap
use rule index as index_re with:
    input: 
        rules.rename_fasta_header_re.output,
    output: 
        temp(OUTPUT_DIR + "/rep_seqs_requant.mmi")
    log: 
        OUTPUT_DIR + "/logs/requant/index.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/index.txt"

use rule dict as dict_re with:
    input: 
        rules.rename_fasta_header_re.output,
    output: 
        temp(OUTPUT_DIR + "/rep_seqs_requant.dict")
    log: 
        OUTPUT_DIR + "/logs/requant/dict.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/dict.txt"

use rule minimap_rep_seqs as minimap_rep_seqs_re with:
    input:
        fq = OUTPUT_DIR + "/f2requant/{barcode}.fastq",
        mmi = rules.index_re.output,
        dict = rules.dict_re.output,
    output: 
        temp(OUTPUT_DIR + "/requant/mapped/{barcode}.bam")
    log: 
        OUTPUT_DIR + "/logs/requant/minimap/{barcode}.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/minimap/{barcode}.txt"

use rule sort as sort_re with:
    input: 
        rules.minimap_rep_seqs_re.output
    output: 
        temp(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam")
    log: 
        OUTPUT_DIR + "/logs/requant/sort/{barcode}.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/sort/{barcode}.txt"

use rule samtools_index as samtools_index_re with:        
    input: 
        rules.sort_re.output
    output: 
        temp(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam.bai")
    log: 
        OUTPUT_DIR + "/logs/requant/index/{barcode}.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/index/{barcode}.txt"

def get_qout_re(wildcards, type_o):
    barcodes = col_info_requant(wildcards, "fq")
    # extract basename & strip extension
    barcodes = [os.path.basename(barcode).split(".")[0] for barcode in barcodes]
    if type_o == "bam":
        output = expand(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam", barcode=barcodes)
    elif type_o == "bai":
        output = expand(OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam.bai", barcode=barcodes)
    elif type_o == "count":
        output = expand(OUTPUT_DIR + "/requant/mapped/{barcode}.count", barcode=barcodes)
    else:
        raise ValueError("type_o must be 'bam', 'bai', or 'count'")
    return output

# biom format header
use rule rowname_kOTU as rowname_kOTU_re with:
    input:
        bam = lambda wildcards: get_qout_re(wildcards, "bam"),
        bai = lambda wildcards: get_qout_re(wildcards, "bai"),
    output: 
        temp(OUTPUT_DIR + "/requant/rowname_seqs")
    log: 
        OUTPUT_DIR + "/logs/requant/rowname_kOTU.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/rowname_kOTU.txt"
       
use rule seqs_count as seqs_count_re with:
    input:
        bam = OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam",
        bai = OUTPUT_DIR + "/requant/mapped/{barcode}.sorted.bam.bai"
    output: 
        temp(OUTPUT_DIR + "/requant/mapped/{barcode}.count")
    log: 
        OUTPUT_DIR + "/logs/requant/seqs_count/{barcode}.log"
    benchmark: 
        OUTPUT_DIR + "/benchmarks/requant/seqs_count/{barcode}.txt"

use rule count_matrix as count_matrix_re with:
    input:
        rowname_seqs = rules.rowname_kOTU_re.output,
        seqs_count = lambda wildcards: get_qout_re(wildcards, "count"),
    output: 
        OUTPUT_DIR + "/count_matrix_requant.tsv"

def col_info_rep(requant = False, chimera_check = True):
    check_val("requant", requant, bool)
    if requant:
        fo = chimeraF(chimera_check)
    else:
        fo = ["count_matrix_requant.tsv", "rep_seqs_requant.fasta"]
    return fo