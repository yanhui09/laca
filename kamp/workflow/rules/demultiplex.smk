rule guppy:
    input: INPUT_DIR
    output: temp(touch(".guppy"))
    log: "logs/demultiplex_guppy.log"
    benchmark: "benchmarks/demultiplex_guppy.txt"
    threads: config["threads"]["large"]
    params:
        barcode_kits=config["barcode_kits"],
        dir=os.path.join(os.getcwd(), "demultiplexed_guppy"),
    shell: "{workflow.basedir}/resources/ont-guppy/bin/guppy_barcoder -i {input} -s {params.dir} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}"

# consider add minibar? https://github.com/calacademy-research/minibar
# demultiplex? https://github.com/jfjlaros/demultiplex
#rule minibar:
#    input: INPUT_DIR
#    output: temp(touch(".minibar"))
#    conda: "../envs/minibar.yaml"
#    log: "logs/demultiplex_minibar.log"
#    benchmark: "benchmarks/demultiplex_minibar.txt"
#    threads: config["threads"]["large"]
#    params:
#        dir = os.path.join(os.getcwd(), "demultiplexed_minibar"),
#        args=config["minibar_args"],
#    shell: 
#        """
#        mkdir -p {params.dir}/unclassified {params.dir}/mult
#        # create barcode list
#        cut -f1 {workflow.basedir}/resources/data/index.txt| sed 1d > {params.dir}/barcodes.txt 2> {log}
#        for i in {input}/*; do 
#            # skip if without .fastq.gz, .fq.gz, .fastq, or .fq
#            if [[ $i != *.fastq.gz && $i != *.fq.gz && $i != *.fastq && $i != *.fq ]]; then
#                continue
#            fi
#            batch_file=$(basename $i)
#            batch_id=${{batch_file%.*}}
#            mkdir {params.dir}/$batch_id
#            echo $batch_id >> {log}
#            python {workflow.basedir}/scripts/minibar.py {workflow.basedir}/resources/data/index.txt $i \
#            -F -T -P {params.dir}/$batch_id/ {params.args} 2>> {log}
#            # sample in barcode list in a dir with batchid
#            while read p; do
#                # if file exists, mkdir and move file
#                if [ -f {params.dir}/$batch_id/$p.fastq ]; then
#                    # if dir not exists, mkdir
#                    if [ ! -d {params.dir}/$p ]; then
#                        mkdir {params.dir}/$p 
#                    fi
#                    mv {params.dir}/$batch_id/$p.fastq {params.dir}/$p/$batch_id.fastq
#                fi
#            done < {params.dir}/barcodes.txt
#            mv {params.dir}/$batch_id/unk.fastq {params.dir}/unclassified/$batch_id.fastq
#            mkdir {params.dir}/mult/$batch_id
#            mv {params.dir}/$batch_id/*.fastq {params.dir}/mult/$batch_id
#            # rm temp dir
#            rmdir {params.dir}/$batch_id
#        done
#        """

# parallel minibar
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
        args = config["minibar_args"],
    log: "logs/minibar/{basecalled_fq}.log"
    benchmark: "benchmarks/minibar/{basecalled_fq}.txt"
    shell: 
        """
        mkdir -p {params.dir}
        python {workflow.basedir}/scripts/minibar.py {workflow.basedir}/resources/data/index.txt {input} \
        -F -T -P {params.dir}/ {params.args} 2> {log}
        """

rule collect_minibar_batch:
    input: expand("demultiplexed_minibar/.demult_{basecalled_fq}", basecalled_fq=get_basecalled_fqs(INPUT_DIR)) 
    output: temp(touch(".minibar"))
    params:
        dir = os.path.join(os.getcwd(), "demultiplexed_minibar"),
    log: "logs/collect_minibar_batch.log"
    benchmark: "benchmarks/collect_minibar_batch.txt"
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
            mv {params.dir}/$batch_id/unk.fastq {params.dir}/unclassified/$batch_id.fastq
            mkdir {params.dir}/mult/$batch_id
            mv {params.dir}/$batch_id/*.fastq {params.dir}/mult/$batch_id
            # rm temp dir
            rmdir {params.dir}/$batch_id
            rm {input} -f
        done
        """

# choose demultiplexer
def get_demult(demult="guppy", dir=False):
    # if demult != "guppy" | "minibar", raise value error
    if demult != "guppy" and demult != "minibar":
        raise ValueError("Demultiplexer not recognized. Choose guppy or minibar in config.")
    demult_dir = os.path.join(os.getcwd(), "demultiplexed_" + demult)
    demult_flag = "." + demult
    if dir:
        return demult_dir
    else:
        return demult_flag

checkpoint demultiplex_check:
    input: ancient(get_demult(demult=config["demultiplex"], dir=False))
    #input: ancient(".minibar")
    output: directory("demultiplexed")
    log: "logs/demultiplex_check.log"
    benchmark: "benchmarks/demultiplex_check.txt"
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

# collect demultiplexed files
rule collect_fastq:
    input:  "demultiplexed/{barcode}"
    output: temp("qc/{barcode}.fastq")
    shell: "cat {input}/*.fastq > {output}"

def get_demultiplexed(wildcards):
    barcodes = glob_wildcards(checkpoints.demultiplex_check.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    return sorted(set(barcodes))
