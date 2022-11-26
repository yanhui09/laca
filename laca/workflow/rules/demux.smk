# None type check in Smakefile, at least one of basecalled_dir and demultiplexed_dir not None
def check_demux_dir(dir_path = config["demultiplexed_dir"]):
    if dir_path is not None:
        if not os.path.isdir(dir_path):
            raise ValueError("\n  'demultiplexed_dir' ({}) in config not found.\n".format(dir_path))
        else:
            if not os.listdir(dir_path):
                raise ValueError("\n  'demultiplexed_dir' ({}) in config is empty.\n".format(dir_path))
            else:
                # shall contain at least one barcode sub-directory, {barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq
                glob_pattern = os.path.join(dir_path, "[a-zA-Z]*[0-9]*", "*.fastq")
                import glob
                if not glob.glob(glob_pattern):
                    raise ValueError("\n  'demultiplexed_dir' ({}) in config does not contain barcode sub-directories.\n  A barcode folder contains unzipped fastq files, named as [a-zA-Z]+[0-9]+.\n".format(dir_path))
check_demux_dir()

rule get_demux_external:
    output: temp(directory("demux_external"))
    params:
        indir = config["demultiplexed_dir"] 
    shell: "cp -r {params.indir} {output}"

def check_basecall_dir(dir_path = config["basecalled_dir"]):
    if dir_path is not None:
        if not os.path.isdir(dir_path):
            raise ValueError("\n  'basecalled_dir' ({}) in config not found.\n".format(dir_path))
        else:
            if not os.listdir(dir_path):
                raise ValueError("\n  'basecalled_dir' ({}) in config is empty.\n".format(dir_path))
check_basecall_dir()

rule guppy:
    # need to bind INPUT_DIR if not in workdir
    input: config["basecalled_dir"]
    output: temp(directory("demux_guppy"))
    singularity: "docker://genomicpariscentre/guppy:3.3.3"
    log: "logs/demultiplex/guppy.log"
    benchmark: "benchmarks/demultiplex/guppy.txt"
    threads: config["threads"]["large"]
    params:
        barcode_kits=config["guppy"]["barcode_kits"],
    shell: 
        """
        guppy_barcoder -i {input} -s {output} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}
        """
    
def get_basecalled_fqs(dir_path=config["basecalled_dir"]):
    basecalled_fqs = []
    suffixes = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
    for fq in os.listdir(dir_path):
        if fq.endswith(suffixes):
            basecalled_fqs.append(fq)
    return basecalled_fqs

rule minibar:
    output: temp(touch("minibar/.demux_{basecalled_fq}"))
    conda: "../envs/minibar.yaml"
    params:
        indir = config["basecalled_dir"], 
        fq = "{basecalled_fq}",
        outdir = lambda wildcards: os.path.join(os.getcwd(), "minibar", os.path.splitext(wildcards.basecalled_fq)[0]),
        args = config["minibar"]["args"],
    log: "logs/demultiplex/minibar/{basecalled_fq}.log"
    benchmark: "benchmarks/demultiplex/minibar/{basecalled_fq}.txt"
    shell: 
        """
        mkdir -p {params.outdir}
        python {workflow.basedir}/scripts/minibar.py {workflow.basedir}/resources/data/index.txt {params.indir}/{params.fq} \
        -F -T -P {params.outdir}/ {params.args} 2> {log}
        """

rule collect_minibars:
    input: expand("minibar/.demux_{basecalled_fq}", basecalled_fq=get_basecalled_fqs())
    output: temp(directory("demux_minibar"))
    params:
        indir = os.path.join(os.getcwd(), "minibar"),
    log: "logs/demultiplex/minibar/collect_batch.log"
    benchmark: "benchmarks/demultiplex/minibar/collect_batch.txt"
    shell:
        """
        mv {params.indir} {output}
        mkdir -p {output}/unclassified {output}/mult
        # create barcode list
        cut -f1 {workflow.basedir}/resources/data/index.txt | sed 1d > {output}/barcodes.txt 2> {log}
        for i in {output}/*/; do 
            if [ "$i" = "{output}/unclassified/" ] || [ "$i" = "{output}/mult/" ]
            then
                continue
            fi
            # sample in barcode list in a dir with batchid
            batch_id=$(basename $i)
            while read p; do
                if [ -f {output}/$batch_id/$p.fastq ]; then
                    # if dir not exists, mkdir and move file
                    if [ ! -d {output}/$p ]; then
                        mkdir {output}/$p 
                    fi
                    mv {output}/$batch_id/$p.fastq {output}/$p/$batch_id.fastq
                fi
            done < {output}/barcodes.txt
            # if file exists, mv to unclassified
            if [ -f {output}/$batch_id/unk.fastq ]; then
               mv {output}/$batch_id/unk.fastq {output}/unclassified/$batch_id.fastq
            fi
            mkdir {output}/mult/$batch_id
            # if .fastq file exists, mv to mult
            if [ -n "$(ls -A {output}/$batch_id/*.fastq 2>/dev/null)" ]; then
               mv {output}/$batch_id/*.fastq {output}/mult/$batch_id
            fi
            # rm temp dir
            rmdir {output}/$batch_id
            rm {input} -f
        done
        """
      
# get demux input
def get_demux(demux=config["demuxer"], demux_external=config["demultiplexed_dir"]):
    if demux_external is not None:
        return rules.get_demux_external.output
    else:
        if demux != "guppy" and demux != "minibar":
            raise ValueError("Demultiplexer not recognized. Choose guppy or minibar in config.")
        return "demux_" + demux

checkpoint demux_check:
    input: get_demux()
    output: directory("demultiplexed")
    log: "logs/demultiplex/check.log"
    benchmark: "benchmarks/demultiplex/check.txt"
    params:
        nreads_m=config["nreads_m"],
    shell: 
        """
        mv {input} {output} > {log} 2>&1
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
    conda: "../envs/seqkit.yaml"
    threads: 1
    log: "logs/demultiplex/collect_fastq/{barcode}.log"
    benchmark: "benchmarks/demultiplex/collect_fastq/{barcode}.txt"
    shell: "cat {input}/*.fastq | seqkit rename -w0 -o {output} 2>> {log}"

def get_demux_barcodes(wildcards):
    barcodes = glob_wildcards(checkpoints.demux_check.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    return sorted(set(barcodes))
