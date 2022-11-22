# None type check in Smakefile, at least one of basecalled_dir and demultiplexed_dir not None
def check_demult_dir(dir_path = config["demultiplexed_dir"]):
    if dir_path is not None:
        if not os.path.isdir(dir_path):
            raise ValueError("\n  Directory {} not found.\n".format(dir_path))
        else:
            if not os.listdir(dir_path):
                raise ValueError("\n  Directory {} is empty.\n".format(dir_path))
            else:
                # shall contain at least one barcode sub-directory, {barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq
                glob_pattern = os.path.join(dir_path, "[a-zA-Z]*[0-9]*", "*.fastq")
                import glob
                if not glob.glob(glob_pattern):
                    raise ValueError("\n  Directory {} does not contain barcode sub-directories.\n  A barcode folder contains unzipped fastq files, named as [a-zA-Z]+[0-9]+.\n".format(dir_path))
check_demult_dir()

rule get_demult_external:
    output: touch(".external_DONE")
    params:
        indir = config["demultiplexed_dir"],
        outdir = os.path.join(os.getcwd(), "demultiplexed_external"),
    shell: "cp -r {params.indir} {params.outdir}"

def check_basecall_dir(dir_path = config["basecalled_dir"]):
    if dir_path is not None:
        if not os.path.isdir(dir_path):
            raise ValueError("\n  Directory {} not found.\n".format(dir_path))
        else:
            if not os.listdir(dir_path):
                raise ValueError("\n  Directory {} is empty.\n".format(dir_path))
check_basecall_dir()

rule guppy:
    # need to bind INPUT_DIR if not in workdir
    output: touch(".guppy_DONE")
    singularity: "docker://genomicpariscentre/guppy:3.3.3"
    log: "logs/demultiplex/guppy.log"
    benchmark: "benchmarks/demultiplex/guppy.txt"
    threads: config["threads"]["large"]
    params:
        barcode_kits=config["guppy"]["barcode_kits"],
        indir=config["basecalled_dir"],
        outdir=os.path.join(os.getcwd(), "demultiplexed_guppy"),
    shell: 
        """
        guppy_barcoder -i {params.indir} -s {params.outdir} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}
        """
    
def get_basecalled_fqs(dir_path=config["basecalled_dir"]):
    basecalled_fqs = []
    suffixes = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
    for fq in os.listdir(dir_path):
        if fq.endswith(suffixes):
            basecalled_fqs.append(fq)
    return basecalled_fqs

rule minibar_batch:
    output: temp(touch("demultiplexed_minibar/.demult_{basecalled_fq}"))
    conda: "../envs/minibar.yaml"
    params:
        indir=config["basecalled_dir"],
        fq = "{basecalled_fq}",
        outdir = lambda wildcards: os.path.join(os.getcwd(), "demultiplexed_minibar", os.path.splitext(wildcards.basecalled_fq)[0]),
        args = config["minibar"]["args"],
    log: "logs/demultiplex/minibar/{basecalled_fq}.log"
    benchmark: "benchmarks/demultiplex/minibar/{basecalled_fq}.txt"
    shell: 
        """
        mkdir -p {params.outdir}
        python {workflow.basedir}/scripts/minibar.py {workflow.basedir}/resources/data/index.txt {params.indir}/{params.fq} \
        -F -T -P {params.outdir}/ {params.args} 2> {log}
        """

rule collect_minibar_batch:
    input: expand("demultiplexed_minibar/.demult_{basecalled_fq}", basecalled_fq=get_basecalled_fqs())
    output: touch(".minibar_DONE")
    params:
        outdir = os.path.join(os.getcwd(), "demultiplexed_minibar"),
    log: "logs/demultiplex/minibar/collect_batch.log"
    benchmark: "benchmarks/demultiplex/minibar/collect_batch.txt"
    shell:
        """
        mkdir -p {params.outdir}/unclassified {params.outdir}/mult
        # create barcode list
        cut -f1 {workflow.basedir}/resources/data/index.txt | sed 1d > {params.outdir}/barcodes.txt 2> {log}
        for i in {params.outdir}/*/; do 
            if [ "$i" = "{params.outdir}/unclassified/" ] || [ "$i" = "{params.outdir}/mult/" ]
            then
                continue
            fi
            # sample in barcode list in a dir with batchid
            batch_id=$(basename $i)
            while read p; do
                if [ -f {params.outdir}/$batch_id/$p.fastq ]; then
                    # if dir not exists, mkdir and move file
                    if [ ! -d {params.outdir}/$p ]; then
                        mkdir {params.outdir}/$p 
                    fi
                    mv {params.outdir}/$batch_id/$p.fastq {params.outdir}/$p/$batch_id.fastq
                fi
            done < {params.outdir}/barcodes.txt
            # if file exists, mv to unclassified
            if [ -f {params.outdir}/$batch_id/unk.fastq ]; then
               mv {params.outdir}/$batch_id/unk.fastq {params.outdir}/unclassified/$batch_id.fastq
            fi
            mkdir {params.outdir}/mult/$batch_id
            # if .fastq file exists, mv to mult
            if [ -n "$(ls -A {params.outdir}/$batch_id/*.fastq 2>/dev/null)" ]; then
               mv {params.outdir}/$batch_id/*.fastq {params.outdir}/mult/$batch_id
            fi
            # rm temp dir
            rmdir {params.outdir}/$batch_id
            rm {input} -f
        done
        """
      
# choose demultiplexer
def get_demult(demult=config["demultiplex"], dir=False):
    if demult != "guppy" and demult != "minibar" and demult != "external":
        raise ValueError("Demultiplexer not recognized. Choose guppy, minibar or external in config.")
    demult_dir = os.path.join(os.getcwd(), "demultiplexed_" + demult)
    demult_flag = "." + demult + "_DONE"
    if dir:
        return demult_dir
    else:
        return demult_flag

checkpoint demultiplex_check:
    input: ancient(get_demult(dir=False))
    output: directory("demultiplexed")
    log: "logs/demultiplex/check.log"
    benchmark: "benchmarks/demultiplex/check.txt"
    params:
        indir=get_demult(dir=True),
        nreads_m=config["nreads_m"],
    shell: 
        """
        mv {params.indir} {output} > {log} 2>&1
        # rm shallow sequencing to avoid bardcode bleeding
        mkdir {output}/suspected -p >> {log} 2>&1
        for i in {output}/*/
        do
            if [ "$i" = "{output}/suspected/" ] || [ "$i" = "{output}/unclassified/" ] || [ "$i" = "{output}/mult/" ]
            then
                continue
            fi
            nlines=$(cat $i/* | wc -l)
            nreads=$((nlines / 4))
            if [ $nreads -lt {params.nreads_m} ]
            then
                mv "$i" {output}/suspected/ >> {log} 2>&1
                echo "$i moved to {output}/suspected/ due to shallow sequencing < {params.nreads_m}" >> {log} 2>&1
            fi
        done
        """

rule collect_fastq:
    input:  "demultiplexed/{barcode}"
    output: temp("qc/{barcode}.fastq")
    shell: "cat {input}/*.fastq > {output}"

def get_demultiplexed(wildcards):
    barcodes = glob_wildcards(checkpoints.demultiplex_check.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    return sorted(set(barcodes))
