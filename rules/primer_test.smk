def check_val(var, val, class_type):
    if not isinstance(val, class_type):
        warns = ('\n\t' + str(var) + ' only accepts ' + str(class_type) + ' values.' +
         '\n\t' + str(val) + ' is used in config.yaml file.')
        raise ValueError(warns)

OUTPUT_DIR = "/mnt/md0/UMI16S/primer_test"
# subsample
subs = False
P = 0.9
N = 1000

# link to download database
dict_db = {
    "greengene": "https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_tax_with_hOTUs_99_reps.fasta",
    "silva": "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
}

# linker and primer info
fprimers = config["fprimer"]
rprimers = config["rprimer"]

# reverse complementation
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
fprimersR = {k: revcomp(v) for k, v in fprimers.items()}
rprimersR = {k: revcomp(v) for k, v in rprimers.items()}

# pattern search for umi using cutadapt
# nanopore possibly sequences either strand
def seqs_join(primer1, primer2):
    joined = '-g ' + primer1 + '...' + revcomp(primer2)
    return joined
def linked_pattern(primers1, primers2):
    primers1_values = list(primers1.values())
    primers2_values = list(primers2.values())
    linked = [seqs_join(primer1, primer2) for primer1 in primers1_values for primer2 in primers2_values]
    return ' '.join(linked)
# generate pattern directory
def pattern_dict(fprimers, rprimers):
    f5_dict = {}
    f5_all = ''
    for i in range(len(fprimers.items())):
        for j in range(len(rprimers.items())):
            f5_pattern1 = linked_pattern(dict(list(fprimers.items())[i:i+1]), dict(list(rprimers.items())[j:j+1]))
            f5_pattern2 = linked_pattern(dict(list(rprimers.items())[j:j+1]), dict(list(fprimers.items())[i:i+1]))
            f5_patterns = f5_pattern1 + ' ' + f5_pattern2
            # update dict
            f5_ij = {str(i+1)+str(j+1): f5_patterns}
            f5_dict.update(f5_ij)
            # include combinations in all
            f5_all += f5_patterns + ' '
    f5_alldict = {'all': f5_all}
    f5_dict.update(f5_alldict)
    return f5_dict

f5_dict = pattern_dict(fprimers, rprimers)

def get_val(key, dict_i):
    return dict_i[key]

rule download_markerDB:
    output: temp(OUTPUT_DIR + "/{db}/rep.fna")
    params:
        _dir = OUTPUT_DIR + "/{db}",
        db_links = lambda wc: get_val(wc.db, dict_db),
    log: OUTPUT_DIR + "/logs/{db}/download_markerDB.log"
    benchmark: OUTPUT_DIR + "/benchmarks/{db}/download_markerDB.txt"
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
    output: OUTPUT_DIR + "/{db}/repDNA.fna"
    conda: "../envs/seqkit.yaml"
    log: OUTPUT_DIR + "/logs/{db}/rna2dna.log"
    benchmark: OUTPUT_DIR + "/benchmarks/{db}/rna2dna.txt"
    shell: "seqkit seq -w 0 --rna2dna {input} > {output} 2> {log}"

rule subsample:
    input: rules.rna2dna.output
    output:
        p = temp(OUTPUT_DIR + "/{db}/subsample/rep_p.fna"),
        n = temp(OUTPUT_DIR + "/{db}/subsample/rep.fna"),
    conda: "../envs/seqkit.yaml"
    params:
        p = P,
        n = N,
    log: OUTPUT_DIR + "/logs/{db}/subsample.log"
    benchmark: OUTPUT_DIR + "/benchmarks/{db}/subsample.txt"
    threads: 1
    shell:
        # I don't know why pipe fails, 
        # portion subsampling (required by seqkit sample by number) as temp file instead.
        """
        seqkit sample -p {params.p} -j {threads} {input} -o {output.p} 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -o {output.n} 2>> {log}
        """

def get_raw(subsample, p, n):
    check_val("subsample", subsample, bool)
    check_val("p[seqkit]", p, float)
    check_val("n[seqkit]", n, int)
    if subsample == True:
        return rules.subsample.output.n
    else:
        return OUTPUT_DIR + "/{db}/repDNA.fna"

# trim primers 
rule trim_primers:
    input: get_raw(subs, P, N)
    output: 
        fna = OUTPUT_DIR + "/{db}/{sets}/primers_trimmed.fna",
        json = OUTPUT_DIR + "/{db}/{sets}/cutadapt.json",
    conda: "../envs/cutadapt.yaml"
    params:
        f = lambda wc: get_val(wc.sets, f5_dict),
    log: OUTPUT_DIR + "/logs/{db}/trim_{sets}.log"
    benchmark: OUTPUT_DIR + "/benchmarks/{db}/trim_{sets}.txt"
    threads: config["threads"]["large"]
    shell:
        """
        cutadapt \
            -j {threads} \
            --discard-untrimmed \
            {params.f} \
            -o {output.fna} \
            {input} \
            --json={output.json} \
            > {log} 2>&1
        """

rule primer_check:
    input: expand(OUTPUT_DIR + "/{db}/{sets}/cutadapt.json", db = [k for k in dict_db.keys()], sets = ["11", "12", "21", "22", "31", "32", "41", "42", "all"])
    output: touch(OUTPUT_DIR + "/.primerCHK_DONE")