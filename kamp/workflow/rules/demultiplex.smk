rule guppy:
    # need to bind INPUT_DIR if not in workdir
    input: INPUT_DIR 
    output: touch(".guppy_DONE")
    singularity: "docker://genomicpariscentre/guppy:3.3.3"
    log: "logs/demultiplex/guppy.log"
    benchmark: "benchmarks/demultiplex/guppy.txt"
    threads: config["threads"]["large"]
    params:
        barcode_kits=config["guppy"]["barcode_kits"],
        dir=os.path.join(os.getcwd(), "demultiplexed_guppy"),
    shell: 
        """
        guppy_barcoder -i {input} -s {params.dir} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}
        """
    
def get_basecalled_fqs(fq_dir):
    basecalled_fqs = []
    suffixes = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
    for fq in os.listdir(fq_dir):
        if fq.endswith(suffixes):
            basecalled_fqs.append(fq)
    return basecalled_fqs

rule minibar_batch:
    input: INPUT_DIR + "/{basecalled_fq}"
    output: temp(touch("demultiplexed_minibar/.demult_{basecalled_fq}"))
    conda: "../envs/minibar.yaml"
    params: 
        dir = lambda wildcards: os.path.join(os.getcwd(), "demultiplexed_minibar", os.path.splitext(wildcards.basecalled_fq)[0]),
        args = config["minibar"]["args"],
    log: "logs/demultiplex/minibar/{basecalled_fq}.log"
    benchmark: "benchmarks/demultiplex/minibar/{basecalled_fq}.txt"
    shell: 
        """
        mkdir -p {params.dir}
        python {workflow.basedir}/scripts/minibar.py {workflow.basedir}/resources/data/index.txt {input} \
        -F -T -P {params.dir}/ {params.args} 2> {log}
        """

rule collect_minibar_batch:
    input: expand("demultiplexed_minibar/.demult_{basecalled_fq}", basecalled_fq=get_basecalled_fqs(INPUT_DIR)) 
    output: touch(".minibar_DONE")
    params:
        dir = os.path.join(os.getcwd(), "demultiplexed_minibar"),
    log: "logs/demultiplex/minibar/collect_batch.log"
    benchmark: "benchmarks/demultiplex/minibar/collect_batch.txt"
    shell:
        """
        mkdir -p {params.dir}/unclassified {params.dir}/mult
        # create barcode list
        cut -f1 {workflow.basedir}/resources/data/index.txt | sed 1d > {params.dir}/barcodes.txt 2> {log}
        for i in {params.dir}/*/; do 
            if [ "$i" = "{params.dir}/unclassified/" ] || [ "$i" = "{params.dir}/mult/" ]
            then
                continue
            fi
            # sample in barcode list in a dir with batchid
            batch_id=$(basename $i)
            while read p; do
                if [ -f {params.dir}/$batch_id/$p.fastq ]; then
                    # if dir not exists, mkdir and move file
                    if [ ! -d {params.dir}/$p ]; then
                        mkdir {params.dir}/$p 
                    fi
                    mv {params.dir}/$batch_id/$p.fastq {params.dir}/$p/$batch_id.fastq
                fi
            done < {params.dir}/barcodes.txt
            # if file exists, mv to unclassified
            if [ -f {params.dir}/$batch_id/unk.fastq ]; then
               mv {params.dir}/$batch_id/unk.fastq {params.dir}/unclassified/$batch_id.fastq
            fi
            mkdir {params.dir}/mult/$batch_id
            # if .fastq file exists, mv to mult
            if [ -n "$(ls -A {params.dir}/$batch_id/*.fastq 2>/dev/null)" ]; then
               mv {params.dir}/$batch_id/*.fastq {params.dir}/mult/$batch_id
            fi
            # rm temp dir
            rmdir {params.dir}/$batch_id
            rm {input} -f
        done
        """

# choose demultiplexer
def get_demult(demult="guppy", dir=False):
    if demult != "guppy" and demult != "minibar":
        raise ValueError("Demultiplexer not recognized. Choose guppy or minibar in config.")
    demult_dir = os.path.join(os.getcwd(), "demultiplexed_" + demult)
    demult_flag = "." + demult + "_DONE"
    if dir:
        return demult_dir
    else:
        return demult_flag

checkpoint demultiplex_check:
    input: ancient(get_demult(demult=config["demultiplex"], dir=False))
    output: directory("demultiplexed")
    log: "logs/demultiplex/check.log"
    benchmark: "benchmarks/demultiplex/check.txt"
    params:
        dir=get_demult(demult=config["demultiplex"], dir=True),
        nreads_m=config["nreads_m"],
    shell: 
        """
        mv {params.dir} {output} > {log} 2>&1
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
