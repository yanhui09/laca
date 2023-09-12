# linker and primer info
flinker = config["flinker"]
fprimers = config["fprimer"]
rlinker = config["rlinker"]
rprimers = config["rprimer"]

# reverse complementation
flinkerR = revcomp(flinker)
rlinkerR = revcomp(rlinker)
fprimersR = {k: revcomp(v) for k, v in fprimers.items()}
rprimersR = {k: revcomp(v) for k, v in rprimers.items()}

# pattern search for umi using cutadapt
# nanopore possibly sequences either strand
def seqs_join_umi(linker, primer, reverse=False):
    joined = '-g ' + linker + '...' + primer
    if reverse:
        joined = '-G ' + primer + '...' + linker
    return joined
def linked_pattern_umi(linker, primers,reverse=False):
    linked = {k: seqs_join_umi(linker, v, reverse) for k, v in primers.items()}        
    return ' '.join(v for v in linked.values())

# forward
f_pattern1 = linked_pattern_umi(flinker, fprimers)
f_pattern2 = linked_pattern_umi(rlinker, rprimers)
f_pattern = f_pattern1 + ' ' + f_pattern2
# reverse
r_pattern1 = linked_pattern_umi(flinkerR, fprimersR, reverse=True)
r_pattern2 = linked_pattern_umi(rlinkerR, rprimersR, reverse=True)
r_pattern = r_pattern1 + ' ' + r_pattern2
#---------
localrules: exclude_shallow_umi, umi_loc, concat_umi, umi_check1, check_umi, umi_check2, umi_loc2, cls_umiCon, split_umibin, fq2fa_umi, collect_umiCon_trimmed

use rule prepare_umapclust_fqs as prepare_umiCon_fqs with:
    input: 
        "qc/qfilt/{barcode}.fastq"
    output: 
        temp("umiCon/fqs/{barcode}_0_0_0.fastq")

use rule prepare_umapclust_fqs as prepare_umiCon_fqs_globalclust with:
    input: 
        rules.fqs_split_meshclust.output.members,
    output: 
        temp("umiCon/fqs/{barcode}_{c1}_{c2}_{c3}.fastq")

def get_fqs4umiCon(wildcards, global_clust = config["globalclust_umi_reads"], expanded = False, ids = False, pool = config["pool"]):
    if expanded == False:
        if global_clust == True:
            return rules.prepare_umiCon_fqs_globalclust.output
        else:
            return rules.prepare_umiCon_fqs.output
    else:
        if global_clust == True:
            bc_clss = glob_wildcards(checkpoints.cls_meshclust.get(**wildcards).output[0] + "/{bc_clss}.csv").bc_clss
        else:
            if pool == True:
               bcs = ["pooled"]
            else:
               bcs = get_qced_barcodes(wildcards)
            bc_clss = [bc +"_0_0_0" for bc in bcs]
        if ids == True:
            return bc_clss
        else: 
            return expand("umiCon/fqs/{bc_clss}.fastq", bc_clss=bc_clss)

# Remove UMI-tagged samples in shallow sequencing
checkpoint exclude_shallow_umi:
    input: lambda wc: get_fqs4umiCon(wc, expanded=True)
    output: temp(directory("umiCon/shallow"))
    params:
        min_reads = 20,
    run:
        import shutil
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            num_reads = sum(1 for line in open(i)) / 4
            if num_reads < params.min_reads:
                shutil.move(i, output[0])

def get_qced_barcodes_umi(wildcards):
    bc_clss = get_fqs4umiCon(wildcards, expanded=True, ids=True)
    bc_clss_shallow = glob_wildcards(checkpoints.exclude_shallow_umi.get(**wildcards).output[0]
     + "/{bc_clss}.fastq").bc_clss
    bc_clss_shallow = sorted(set(bc_clss_shallow))
    bc_clss_pass = [b for b in bc_clss if b not in bc_clss_shallow]
    return bc_clss_pass

# trim umi region
rule umi_loc:
    input: lambda wc: get_fqs4umiCon(wc)
    output:
        start = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/start.fastq"),
        end = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/end.fastq"),
    params:
        umi_loc=config["umi"]["loc"]
    shell:
        """
        seqkit subseq -r 1:{params.umi_loc} {input} --quiet 1> {output.start}
        seqkit subseq -r -{params.umi_loc}:-1 {input} --quiet 1> {output.end}
        """

# extract UMI sequences
rule extract_umi:
    input:
        start = rules.umi_loc.output.start,
        end = rules.umi_loc.output.end,
    output:
        umi1 = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi1.fastq"),
        umi2 = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi2.fastq"),
    conda: "../envs/cutadapt.yaml"
    params:
        f = f_pattern,
        r = r_pattern,
        max_err = config["umi"]["cutadapt"]["max_errors"],
        min_overlap = config["umi"]["cutadapt"]["min_overlap"],
        min_len = config["umi"]["len"],
        max_len = config["umi"]["len"],
    log: "logs/umiCon/umiExtract/extract_umi/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/umiCon/umiExtract/extract_umi/{barcode}_{c1}_{c2}_{c3}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        "cutadapt "
        "-j {threads} -e {params.max_err} -O {params.min_overlap} "
        "-m {params.min_len} -M {params.max_len} "
        "--discard-untrimmed "
        "{params.f} "
        "{params.r} "
        "-o {output.umi1} -p {output.umi2} "
        "{input.start} {input.end} > {log} 2>&1"

# combine UMI sequences
rule concat_umi:
    input:
        umi1 = rules.extract_umi.output.umi1,
        umi2 = rules.extract_umi.output.umi2,
    output: temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi.fasta")
    shell: "seqkit concat {input.umi1} {input.umi2} --quiet | seqkit fq2fa -o {output} --quiet"

# extend list with umi_loc fqs
def get_umifile1(wildcards, extend=False):
    bc_clss = get_qced_barcodes_umi(wildcards)
    fs = expand("umiCon/umiExtract/{bc_clss}/umi.fasta", bc_clss=bc_clss)
    if extend == True:
        fs.extend(expand("umiCon/umiExtract/{bc_clss}/{loc}.fastq", bc_clss=bc_clss, loc=["start", "end"]))
        fs.extend(expand("umiCon/fqs/{bc_clss}.fastq", bc_clss=bc_clss))
    return fs     

checkpoint umi_check1:
    input: 
        lambda wc: get_umifile1(wc, extend=True),
        "umiCon/shallow",
    output: temp(directory("umiCon/umiExtract/check1"))
    run:
        import shutil
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        # {input} couldn't be evaluted correctly in `run`, re-evaluating it as walkround solution
        # possibly related to https://github.com/snakemake/snakemake/issues/55
        fs = get_umifile1(wildcards)
        for i in list(fs):
            if os.stat(i).st_size > 0:
                bc_clss = i.split("/")[-2]
                shutil.move(i, output[0] + "/{bc_clss}.fasta".format(bc_clss=bc_clss))

# check UMI pattern
# generate grep pattern for UMI regions
def get_umi_pattern(umi_pattern, length):
    # build a hash table for IUPAC codes
    iupac = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT"
    }
    # rm whitespaces
    umi_pattern = umi_pattern.replace(" ", "")

    # check the length of umi pattern
    n = len(umi_pattern)
    if n != 2*length:
        raise ValueError("UMI pattern length is not equal to your input.\t\nPlease check your `len` and `umi_pattern` in config.yaml .")
    
    grep_pattern = ""
    iloc = 1
    count = 1
    while iloc < n:
        if umi_pattern[iloc] != umi_pattern[iloc-1]:
            # check if the IUPAC code
            if umi_pattern[iloc] not in iupac or umi_pattern[iloc-1] not in iupac:
                raise ValueError("UMI pattern contains non-IUPAC code.\t\nPlease check your `umi_pattern` in config.yaml .")
            grep_pattern += "[" + iupac[umi_pattern[iloc-1]] + "]{" + str(count) + "}" 
            count = 1
        else:
            count += 1

        if iloc == n-1:
            grep_pattern += "[" + iupac[umi_pattern[iloc]] + "]{" + str(count) + "}"
        iloc += 1
    return grep_pattern

rule check_umi:
    input: "umiCon/umiExtract/check1/{barcode}_{c1}_{c2}_{c3}.fasta"
    output: temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12f.fasta")
    params:
        pattern = lambda wc: get_umi_pattern(config["umi"]["pattern"], config["umi"]["len"])
    shell: "grep -B1 -E '{params.pattern}' {input} | sed '/^--$/d' > {output}"

# cluster UMIs
rule cluster_umi:
    input: rules.check_umi.output
    output:
        umi12u = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12u.fasta"),
        umi12cs = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12cs.fasta"),
        umi12uc = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12uc.txt"),
        umi12c = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12c.fasta"),
    conda: "../envs/vsearch.yaml"
    params:
        cl_identity = config["umi"]["cl_identity"],
    log: "logs/umiCon/umiExtract/cluster_umi/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/umiCon/umiExtract/cluster_umi/{barcode}_{c1}_{c2}_{c3}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["default"],
    shell:
        """
        vsearch --fastx_uniques {input} --fastaout {output.umi12u} \
         --sizeout --minuniquesize 1 --relabel umi --strand both > {log} 2>&1

        vsearch --cluster_fast {output.umi12u} --id {params.cl_identity} \
        --centroids {output.umi12cs} --uc {output.umi12uc} --sizein --sizeout \
        --strand both --threads {threads} >> {log} 2>&1

        vsearch --sortbysize {output.umi12cs} --minsize 1 --output {output.umi12c} >> {log} 2>&1
        """

# align the trimmed UMI reads to the potential UMIs, which determine the final UMI clusters
def trim_primers(primers, reverse=False):
    if isinstance(primers, str):
        primers = dict(seqs = primers)
    patterns = {k: '-g ' + v for k, v in primers.items()}
    if reverse:
        patterns = {k: '-G ' + v for k, v in primers.items()}
    return ' '.join(v for v in patterns.values())

fprimer1_trim = trim_primers(flinker)
fprimer2_trim = trim_primers(rlinker)
fprimers_trim = fprimer1_trim + ' ' + fprimer2_trim

rprimer1_trim = trim_primers(fprimersR, reverse=True)
rprimer2_trim = trim_primers(rprimersR, reverse=True)
rprimers_trim = rprimer1_trim + ' ' + rprimer2_trim

rule extract_umip:
    input:
        start = rules.umi_loc.output.start,
        end = rules.umi_loc.output.end,
    output:
        umi1 = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi1p.fastq"),
        umi2 = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi2p.fastq"),
    conda: "../envs/cutadapt.yaml"
    params:
        f = fprimers_trim,
        r = rprimers_trim,
        max_err = config["umi"]["cutadapt"]["max_errors"],
        min_overlap = config["umi"]["cutadapt"]["min_overlap"],
        min_len = config["umi"]["len"],
        max_len = config["umi"]["len"],
    log: "logs/umiCon/umiExtract/extract_umip/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/umiCon/umiExtract/extract_umip/{barcode}_{c1}_{c2}_{c3}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        "cutadapt "
        "-j {threads} -e {params.max_err} -O {params.min_overlap} "
        "-m {params.min_len} -l {params.max_len} "
        "--discard-untrimmed "
        "{params.f} "
        "{params.r} "
        "-o {output.umi1} -p {output.umi2} "
        "{input.start} {input.end} > {log} 2>&1"

use rule concat_umi as concat_umip with:
    input:
        umi1=rules.extract_umip.output.umi1,
        umi2=rules.extract_umip.output.umi2,
    output: 
        temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12p.fasta")

# use bwa to map the potential UMI reads to the UMI ref
rule bwa_umi:
    input: 
        ref = rules.cluster_umi.output.umi12c,
        umip = rules.concat_umip.output,
    output:
        temp(multiext("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12c.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")),
        sai = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12p.sai"),
        sam = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12p.sam"),
    conda: "../envs/bwa.yaml"
    params:
        n = 6,
        F = 4,
    log: "logs/umiCon/umiExtract/bwa_umi/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/umiCon/umiExtract/bwa_umi/{barcode}_{c1}_{c2}_{c3}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: 
        """
        bwa index {input.ref} 2> {log}
        bwa aln -n {params.n} -t {threads} -N {input.ref} {input.umip} > {output.sai} 2>> {log}
        bwa samse -n 10000000 {input.ref} {output.sai} {input.umip} 2>> {log} | \
        samtools view -F {params.F} - > {output.sam} 2>> {log}
        """

rule umi_filter:
    input:
        sam = rules.bwa_umi.output.sam,
        ref = rules.cluster_umi.output.umi12c,
    output: temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi12cf.fasta")
    conda: "../envs/umi.yaml"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        awk \
        -v UMS={input.sam} \
        -v UC={input.ref} \
        '
        (FILENAME == UMS){{
            CLUSTER[$3]++
        }}
        (FILENAME == UC && FNR%2==1){{
          NAME=substr($1,2)
          if (NAME in CLUSTER){{
            if (CLUSTER[NAME]+0 > 2){{
              SIZE=CLUSTER[NAME]
              gsub(";.*", "", NAME)
              print ">" NAME ";size=" SIZE ";"
              getline; print
            }}
          }}
        }}
        ' \
        {input.sam} \
        {input.ref} \
        > {output} 
        """

# extend list with umi_loc fqs
def get_umifile2(wildcards, extend = False):
    bc_clss = glob_wildcards(checkpoints.umi_check1.get(**wildcards).output[0] + "/{bc_clss}.fasta").bc_clss
    fs = expand("umiCon/umiExtract/{bc_clss}/umi12cf.fasta", bc_clss=bc_clss)
    if extend == True:
        fs.extend(expand("umiCon/umiExtract/{bc_clss}/{loc}.fastq", bc_clss=bc_clss, loc=["start", "end"]))
        fs.extend(expand("umiCon/fqs/{bc_clss}.fastq", bc_clss=bc_clss))
    return fs

checkpoint umi_check2:
    input: 
        "umiCon/umiExtract/check1",
        lambda wc: get_umifile2(wc, extend = True),
    output: temp(directory("umiCon/umiExtract/check2"))
    run:
        import shutil
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        fs = get_umifile2(wildcards)
        for i in list(fs):
            if os.stat(i).st_size > 0:
                bc_cls = i.split("/")[-2]
                shutil.move(i, output[0] + "/{bc_cls}.fasta".format(bc_cls=bc_cls))

# rm potential chimera
rule rm_chimera:
    input: "umiCon/umiExtract/check2/{barcode}_{c1}_{c2}_{c3}.fasta",
    output:
        ref = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi_ref.txt"),
        fasta = temp("umiCon/umiExtract/{barcode}_{c1}_{c2}_{c3}/umi_ref.fasta"),
    conda: "../envs/umi.yaml"
    params:
        umip_len = config["umi"]["len"],
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        paste <(cat {input} | paste - - ) \
          <(awk '!/^>/{{print}}' {input} | rev | tr ATCG TAGC) |\
          awk 'NR==FNR {{
              #Format columns
              split($1, a, /[>;]/);
              sub("size=", "", a[3]);
              # Extract UMI1 and UMI2 in both orientations
              s1 = substr($2, 1, {params.umip_len});
              s2 = substr($2, {params.umip_len}+1, 2*{params.umip_len});
              s1rc= substr($3, 1, {params.umip_len});
              s2rc= substr($3, {params.umip_len}+1, 2*{params.umip_len});
              # Register UMI1 size if larger than current highest or if empty
              if ((g1n[s1]+0) <= (a[3]+0) || g1n[s1] == ""){{
                g1n[s1] = a[3];
                g1[s1] = a[2];
              }}
              # Register UMI2 size if larger than current highest or if empty
              if ((g2n[s2]+0) <= (a[3]+0) || g2n[s2] == ""){{
                g2n[s2] = a[3];
                g2[s2] = a[2];
              }}
              # Register UMI1rc size if larger than current highest or if empty
              if ((g1n[s1rc]+0) <= (a[3]+0) || g1n[s1rc] == ""){{
                g1n[s1rc] = a[3];
                g1[s1rc] = a[2];
              }}
              # Register UMI2rc size if larger than current highest or if empty
              if ((g2n[s2rc]+0) <= (a[3]+0) || g2n[s2rc] == ""){{
                g2n[s2rc] = a[3];
                g2[s2rc] = a[2];
              }}
              # Register UMI1 and UMI matches for current UMI
              u[a[2]] = a[3];
              s1a[a[2]] = s1;
              s2a[a[2]] = s2;
              s1arc[a[2]] = s1rc;
              s2arc[a[2]] = s2rc;
            }} END {{
              for (i in u){{
                keep="no";
                if (g1[s1a[i]] == i && g2[s2a[i]] == i && g1[s1arc[i]] == i && g2[s2arc[i]] == i && s1a[i] != s1arc[i]){{
                  keep="yes";
                  print ">"i";"u[i]"\\n"s1a[i]s2a[i] > "{output.fasta}";
                }} else if (s1a[i] == s1arc[i]){{
                  keep="tandem"
                  print ">"i";"u[i]"\\n"s1a[i]s2a[i] > "{output.fasta}";
                }}
                print i, n[i], s1a[i], s2a[i], keep, g1[s1a[i]]"/"g2[s2a[i]]"/"g1[s1arc[i]]"/"g2[s2arc[i]], u[i]
              }}  
            }}' > {output.ref}
        """

# UMI binning
# extract strict UMI region
rule umi_loc2:
    input:
        start = rules.umi_loc.output.start,
        end = rules.umi_loc.output.end,
    output:
        start = temp("umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/reads_umi1.fasta"),
        end = temp("umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/reads_umi2.fasta"),
    params:
        s = config["umi"]["s"],
        e = config["umi"]["e"],
    shell:
        """
        seqkit subseq -r 1:{params.s} {input.start} --quiet | seqkit fq2fa -w0 -o {output.start} --quiet
        seqkit subseq -r -{params.e}:-1 {input.end} --quiet | seqkit fq2fa -w0 -o {output.end} --quiet
        """

# split UMI ref to capture UMI in both orientations
rule split_umi_ref:
    input: rules.rm_chimera.output.fasta
    output:
        umi1 = temp("umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/barcodes_umi1.fasta"),
        umi2 = temp("umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/barcodes_umi2.fasta"),
    conda: "../envs/umi.yaml"
    params:
        umi_len = config["umi"]["len"],
    shell:
        """
        cat {input} <(seqtk seq -r {input} |\
          awk 'NR%2==1{{print $0 "_rc"; getline; print}};') |\
          awk 'NR%2==1{{
               print $0 > "{output.umi1}";
               print $0 > "{output.umi2}";  
             }}
             NR%2==0{{
               print substr($0, 1, {params.umi_len}) > "{output.umi1}";
               print substr($0, length($0) - {params.umi_len} + 1, {params.umi_len})  > "{output.umi2}";  
             }}'
        """

# Map UMI barcode to regions containing UMI
## Important settings:
## -N : diasble iterative search. All possible hits are found.
## -F 20 : Removes unmapped and reverse read matches. Keeps UMIs
##         in correct orientations.
use rule bwa_umi as bwa_umi2 with:
    input: 
        ref = "umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/reads_{umi}.fasta",
        umip = "umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/barcodes_{umi}.fasta",
    output:
        temp(multiext("umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/reads_{umi}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")),
        sai = temp("umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/{umi}.sai"),
        sam = temp("umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/{umi}.sam"),
    params:
        n = 3,
        F = 20,
    log:
        "logs/umiCon/umiBin/bwa_umi2/{barcode}_{c1}_{c2}_{c3}_{umi}.log"
    benchmark:
        "benchmarks/umiCon/umiBin/bwa_umi2/{barcode}_{c1}_{c2}_{c3}_{umi}.txt"
 
rule umi_binning:
    input:
        umi1 = "umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/umi1.sam",
        umi2 = "umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/umi2.sam",
    output:
        bins = "umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/umi_bins.txt",
        stats = "umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}/umi_stats.txt",
    conda: "../envs/umi.yaml"
    params:
        u = config["umi"]["u"],
        U = config["umi"]["U"],
        O = config["umi"]["O"],
        N = config["umi"]["N"],
        S = config["umi"]["S"],
    log: "logs/umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}.log"
    benchmark: "benchmarks/umiCon/umiBin/{barcode}_{c1}_{c2}_{c3}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        awk \
        -v UME_MATCH_ERROR="{params.u}" \
        -v UME_MATCH_ERROR_SD="{params.U}" \
        -v RO_FRAC="{params.O}" \
        -v MAX_BIN_SIZE="{params.N}"  \
        -v BIN_CLUSTER_RATIO="{params.S}" \
        '
        NR==1 {{
          print "[" strftime("%T") "] ### Read-UMI match filtering ###" > "/dev/stderr";
          print "[" strftime("%T") "] Reading UMI1 match file..." > "/dev/stderr";
        }}
        # Read UMI match file
        NR==FNR {{
          # Extract data from optional fields
          for (i=12; i <= NF; i++){{
            # Find NM field and remove prefix (primary hit err)
            if($i ~ /^NM:i:/){{sub("NM:i:", "", $i); perr = $i}};
            # Find secondary hit field, remove prefix and split hits
            if($i ~ /^XA:Z:/){{sub("XA:Z:", "", $i); split($i, shits, ";")}};
          }}
          # Add primary mapping to hit list
          err1[$1][$3]=perr;
          # Add secondary mapping to hit list
          #Iterate over each hit
          for (i in shits){{
            # Split hit in subdata (read, pos, cigar, err)  
            split(shits[i], tmp, ",");
            # Add hit if non-empty, not seen before and not target reverse strand
            if (tmp[1] != "" && !(tmp[1] in err1[$1]) && tmp[2] ~ "+"){{
              err1[$1][tmp[1]] = tmp[4];
            }}
          }}
          next;
        }}
        FNR==1 {{
         print "[" strftime("%T") "] Reading UMI2 match file..." > "/dev/stderr";
        }}
        # Read UMI match file
        {{
          # Extract data from optional fields
          for (i=12; i <= NF; i++){{
            # Find NM field and remove prefix (primary hit err)
            if($i ~ /^NM:i:/){{sub("NM:i:", "", $i); perr = $i}};
            # Find secondary hit field and remove prefix
            if($i ~ /^XA:Z:/){{sub("XA:Z:", "", $i); split($i, shits, ";")}};
          }}
          # Add primary mapping to hit list
          err2[$1][$3]=perr;
          # Add secondary mapping to hit list
          # Split list of hits 
          #Iterate over each hit
          for (i in shits){{
            # Split hit in subdata (read, pos, cigar, err)
            split(shits[i], tmp, ",");
            # Add hit if non-empty, not seen before and not target reverse strand
            if (tmp[1] != "" && !(tmp[1] in err2[$1]) && tmp[2] ~ "+"){{
              err2[$1][tmp[1]] = tmp[4];
            }}
          }}
        }} END {{
          print "[" strftime("%T") "] UMI match filtering..." > "/dev/stderr"; 
          # Filter reads based on UMI match error
          for (umi in err1){{    
            for (read in err1[umi]){{
              # Define vars
              e1 = err1[umi][read];
              e2 = err2[umi][read];
              # Filter reads not matching both UMIs
              if (e1 != "" && e2 != ""){{
                # Filter based on mapping error 
                if (e1 + e2 <= 6 && e1 <= 3 && e2 <= 3){{
                  # Add read to bin list or replace bin assignment if error is lower
                  if (!(read in match_err)){{
                    match_umi[read] = umi;
                    match_err[read] = e1 + e2;
                  }} else if (match_err[read] > e1 + e2 ){{
                    match_umi[read] = umi;
                    match_err[read] = e1 + e2;
                  }} 
                }}
              }}
            }}
          }}
          print "[" strftime("%T") "] Read orientation filtering..." > "/dev/stderr";
          # Count +/- strand reads
          for (s in match_umi){{
            UM=match_umi[s]
            sub("_rc", "", UM)
            # Read orientation stats
            ROC=match(match_umi[s], /_rc/)
            if (ROC != 0){{
              umi_ro_plus[UM]++
              roc[s]="+"
            }} else {{
              umi_ro_neg[UM]++
              roc[s]="-"
            }}
            # Count reads per UMI bin
            umi_n_raw[UM]++;
          }}
          
          # Calculate read orientation fraction
          for (u in umi_ro_plus){{
            # Check read orientation fraction
            if (umi_ro_plus[u] > 1 && umi_ro_neg[u] > 1){{
              if (umi_ro_plus[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){{
                rof_check[u]="rof_subset"
                rof_sub_neg_n[u] = umi_ro_plus[u]*(1/RO_FRAC-1)
                rof_sub_pos_n[u] = rof_sub_neg_n[u]
              }} else if (umi_ro_neg[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){{
                rof_check[u]="rof_subset"
                rof_sub_neg_n[u]=umi_ro_neg[u]*(1/RO_FRAC-1)
                rof_sub_pos_n[u]=rof_sub_neg_n[u]
              }} else {{
                rof_check[u]="rof_ok"
                rof_sub_neg_n[u]=MAX_BIN_SIZE
                rof_sub_pos_n[u]=MAX_BIN_SIZE
              }}
            }} else {{
              rof_check[u]="rof_fail"
            }}
          }}
          
          # Subset reads
          for (s in match_umi){{
            UMI_NAME=match_umi[s]
            sub("_rc", "", UMI_NAME)
            if(roc[s] == "+"){{
              if(rof_sub_pos_n[UMI_NAME]-- > 0){{
                ror_filt[s]=UMI_NAME
              }}
            }} else if (roc[s] == "-"){{
              if(rof_sub_neg_n[UMI_NAME]-- > 0){{
                ror_filt[s]=UMI_NAME
              }}
            }}
          }}
          print "[" strftime("%T") "] UMI match error filtering..." > "/dev/stderr";
          # Calculate UME stats
          for (s in ror_filt){{
            UM=ror_filt[s]
            # Count matching reads
            umi_n[UM]++;
            # UMI match error stats
            umi_me_sum[UM] += match_err[s]
            umi_me_sq[UM] += (match_err[s])^2
          }}
          # Check UMI match error
          for (u in umi_n){{
            UME_MEAN[u] = umi_me_sum[u]/umi_n[u]
            UME_SD[u] = sqrt((umi_me_sq[u]-umi_me_sum[u]^2/umi_n[u])/umi_n[u])
            if (UME_MEAN[u] > UME_MATCH_ERROR || UME_SD[u] > UME_MATCH_ERROR_SD){{
              ume_check[u] = "ume_fail"
            }} else {{
              ume_check[u] = "ume_ok"
            }}
          }}
          print "[" strftime("%T") "] UMI bin/cluster size ratio filtering..." > "/dev/stderr";
          for (u in umi_n){{
            CLUSTER_SIZE=u
            sub(".*;", "", CLUSTER_SIZE)
            bcr[u]=umi_n_raw[u]/CLUSTER_SIZE
            if (bcr[u] > BIN_CLUSTER_RATIO){{
              bcr_check[u] = "bcr_fail"
            }} else {{
              bcr_check[u] = "bcr_ok"
            }}
          }}
          # Print filtering stats
          print "umi_name", "read_n_raw", "read_n_filt", "read_n_plus", "read_n_neg", \
            "read_max_plus", "read_max_neg", "read_orientation_ratio", "ror_filter", \
            "umi_match_error_mean", "umi_match_error_sd", "ume_filter", "bin_cluster_ratio", \
            "bcr_filter" > "{output.stats}"
          for (u in umi_n){{
            print u, umi_n_raw[u], umi_n[u], umi_ro_plus[u], umi_ro_neg[u], \
              rof_sub_pos_n[u] + umi_ro_plus[u], rof_sub_neg_n[u] + umi_ro_neg[u], rof_check[u], \
              UME_MEAN[u], UME_SD[u], ume_check[u], bcr[u], bcr_check[u]\
              > "{output.stats}"
          }}
        
          print "[" strftime("%T") "] Print UMI matches..." > "/dev/stderr"; 
          for (s in ror_filt){{
            UMI_NAME=ror_filt[s]
            if( \
                ume_check[UMI_NAME] == "ume_ok" && \
                rof_check[UMI_NAME] == "rof_ok" && \
                bcr_check[UMI_NAME] == "bcr_ok" \
            ){{print UMI_NAME, s, match_err[s]}}
          }}
          # Print to terminal
          print "[" strftime("%T") "] Done." > "/dev/stderr"; 
        }}
        ' {input.umi1} {input.umi2} > {output.bins} 2> {log}
        """

def get_filt_umi(wildcards):
    bc_clss = glob_wildcards(checkpoints.umi_check2.get(**wildcards).output[0] + "/{bc_clss}.fasta").bc_clss
    return expand("umiCon/fqs/{bc_clss}.fastq", bc_clss=bc_clss)

# split reads by umi bins
checkpoint cls_umiCon:
    input: 
        "umiCon/umiExtract/check2",
        lambda wc: get_filt_umi(wc),
        clss = lambda wc: expand("umiCon/umiBin/{bc_clss}/umi_bins.txt", bc_clss=glob_wildcards(checkpoints.umi_check2.get(**wc).output[0] + "/{bc_clss}.fasta").bc_clss),
    output: directory("umiCon/clusters")
    params:
        min_size = config["min_support_reads"],
    run:
        import pandas as pd
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input.clss):
            bc_cls = i.split("/")[-2]
            df = pd.read_csv(i, sep=" ", header=None, usecols=[0,1], names=["umi","read"])
            df.umi = [x.split(";")[0] for x in df.umi]
            for umi in set(df.umi):
                umi_reads = df.loc[df.umi == umi, "read"]
                if len(umi_reads) >= params.min_size:
                    umi_reads.to_csv(output[0] + "/" + bc_cls + str(umi).replace("umi","cand") + ".csv", header=False, index=False)
        
use rule fqs_split_isONclust as split_umibin with:
    input: 
        cluster = "umiCon/clusters/{barcode}_{c1}_{c2}_{c3}cand{cand}.csv",
        fqs = lambda wc: get_fqs4umiCon(wc),
    output:
        temp("umiCon/split/{barcode}_{c1}_{c2}_{c3}cand{cand}.fastq"),

# find the seed read
rule fq2fa_umi:
    input: rules.split_umibin.output
    output: temp("umiCon/polish/{barcode}_{c1}_{c2}_{c3}cand{cand}/split.fna")
    shell: "seqkit fq2fa {input} -o {output} --quiet"

rule get_seedfa:
    input: rules.fq2fa_umi.output
    output: 
        centroids = temp("umiCon/polish/{barcode}_{c1}_{c2}_{c3}cand{cand}/centroids.fna"),
        seed = temp("umiCon/polish/{barcode}_{c1}_{c2}_{c3}cand{cand}/minimap2/raw.fna")
    conda: "../envs/vsearch.yaml"
    log: "logs/umiCon/{barcode}_{c1}_{c2}_{c3}cand{cand}/get_centroids.log"
    benchmark: "benchmarks/umiCon/{barcode}_{c1}_{c2}_{c3}cand{cand}/get_centroids.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        vsearch --cluster_fast {input} -id 0.75 --sizeout --strand both \
        --centroids {output.centroids} --threads {threads} > {log} 2>&1

        vsearch --sortbysize {output.centroids} --topn 1 --relabel seed \
        --output {output.seed} --threads {threads} >> {log} 2>&1
        """

# polish follows rules in consensus.smk

# trim umi adaptors 
rule trim_adaptors_umi:
    input: "umiCon/umiCon.fna"
    output: 
        trimmed = temp("umiCon/umiCon_trimmedF.fna"),
        untrimmed = temp("umiCon/umiCon_untrimmedF.fna"),
    conda: "../envs/cutadapt.yaml"
    params:
        f = f5_pattern1,
        e = 0.1,
        O = 3,
        m = config["umi"]["seqkit"]["min_len"],
        M = config["umi"]["seqkit"]["max_len"],
        action = "retain",
    threads: config["threads"]["normal"]
    log: "logs/umiCon/trim_adaptors_umi.log"
    benchmark: "benchmarks/umiCon/trim_adaptors_umi.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        cutadapt \
        --action={params.action} \
        -j {threads} \
        -e {params.e} -O {params.O} -m {params.m} -M {params.M} \
        {params.f} \
        --untrimmed-output {output.untrimmed} \
        -o {output.trimmed} \
        {input} \
        > {log} 2>&1
        """

use rule trim_adaptors_umi as trim_adaptors_umiR with:
    input: 
        rules.trim_adaptors_umi.output.untrimmed
    output: 
        trimmed = temp("umiCon/umiCon_trimmedR.fna"),
        untrimmed = "umiCon/umiCon_untrimmed.fna"
    params:
        f = f5_pattern2,
        e = 0.1,
        O = 3,
        m = config["umi"]["seqkit"]["min_len"],
        M = config["umi"]["seqkit"]["max_len"],
        action = "retain",
    log: 
        "logs/umiCon/trim_primers_umiR.log"
    benchmark: 
        "benchmarks/umiCon/trim_primers_umiR.txt"

# reverse complement for reverse strand
use rule revcomp_fq as revcomp_fq_umi with:
    input: 
        rules.trim_adaptors_umiR.output.trimmed
    output: 
        temp("umiCon/umiCon_trimmedR_revcomp.fna")

rule collect_umiCon_trimmed:
    input: 
        rules.trim_adaptors_umi.output.trimmed,
        rules.revcomp_fq_umi.output
    output: "umiCon/umiCon_trimmed.fna"
    run: 
        with open(output[0], "w") as out:
            for i in input:
                with open(i, "r") as inp:
                    for line in inp:
                        out.write(line)