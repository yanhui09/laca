from .log import logger

import os
from ruamel.yaml import YAML

def init_conf(
    bascdir,
    demuxdir,
    merge,
    merge_parent,
    dbdir,
    workdir,
    config="config.yaml",
    demuxer = "guppy",
    nreads_m = 1000,
    no_pool=False,
    subsample=False,
    no_chimera_filt=False,
    no_primer_check=False,
    cluster=["isONclust", "meshclust"],
    consensus=["kmerCon"],
    quant=["seqid"],
    uchime=False,
    jobs_m=2,
    jobs_M=6,
    ont=False,
    isoseq=False,
    longumi=False,
    simulate=False,
):
    """
    Reads template config file with comments from ./template_config.yaml
    updates it by the parameters provided.
    Args:
        bascdir (str): path to a directory of basecalled fastq files
        demuxdir (str): path to a directory of demultiplexed fastq files
        merge (list): list of the working directories of LACA runs to merge
        merge_parent (str): path to the parent of the working directories of LACA runs to merge
        dbdir (str): path to the taxonomy database
        workdir (str): path to the working directory
        config (str): the config filename
        demuxer (str): the demultiplexer [default: "guppy"]
        nreads_m (int): minimum number of reads for the demultiplexed fastqs
        no_pool (bool): if True, do not pool the reads [default: False]
        subsample (bool): if True, subsample the reads [default: False]
        no_chimera_filt (bool): if True, do not filter possible chimeric reads [default: False]
        no_primer_check (bool): if True, do not check primer pattern [default: False]
        cluster (list): list of methods to cluster reads (isONclust, umapclust, meshclust) [default: ["isONclust", "meshclust"]]
        consensus (list): list of methods to generate consensus (kmerCon, miniCon, isoCon, umiCon) [default: ["kmerCon"]]
        quant (list): list of methods to create abundance matrix (seqid, minimap2) [default: ["seqid"]]
        uchime (bool): if True, filter possible chimeric OTUs by vsearch [default: False]
        jobs_m (int): number of jobs for common tasks [default: 2]
        jobs_M (int): number of jobs for threads-dependent tasks [default: 6]
        ont (bool): if True, use template for ONT reads [default: False]
        isoseq (bool): if True, use template for PacBio isoseq reads [default: False]
        longumi (bool): if True, use primer design from longumi paper (https://doi.org/10.1038/s41592-020-01041-y) [default: False]       
    """
    os.makedirs(dbdir, exist_ok=True)
    os.makedirs(workdir, exist_ok=True)
    
    yaml = YAML()
    template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workflow/config_template.yaml")
    
    with open(template_conf_file) as template_conf:
        conf = yaml.load(template_conf)
    
    # no-flags can be overwritten    
    conf["pool"] = not no_pool
    conf["chimera_filt"] = not no_chimera_filt
    conf["primer_check"] = not no_primer_check
    
    if ont == True:
        # yacrd
        conf["yacrd"]["minimap2"]["g"] = 500
        conf["yacrd"]["c"] = 4
        conf["yacrd"]["n"] = 0.4
        # isONclustCon
        conf["isONclust"]["k"] = 13
        conf["isONclust"]["w"] = 20
        # isoCon with isONcorrect
        conf["isoCon"]["run_isONcorrect"] = True
        # minimap2
        conf["minimap2"]["x_ava"] = "ava-ont"
        conf["minimap2"]["x_map"] = "map-ont"
        # racon medaka
        conf["racon"]["iter"] = 2
        conf["medaka"]["iter"] = 1
        # UMI
        conf["umi"]["s"] = 90
        conf["umi"]["e"] = 90
        conf["umi"]["u"] = 3.5
        conf["umi"]["O"] = 0.2
        # simulate
        conf["simulate"]["badread"]["error_model"] = "nanopore2020"
        conf["simulate"]["badread"]["qscore_model"] = "nanopore2020"
    
    if isoseq == True:
        # yacrd
        conf["yacrd"]["minimap2"]["g"] = 5000
        conf["yacrd"]["c"] = 3
        conf["yacrd"]["n"] = 0.4
        # isONclustCon
        conf["isONclust"]["k"] = 15
        conf["isONclust"]["w"] = 50
        # isoCon without isONcorrect
        conf["isoCon"]["run_isONcorrect"] = False
        # minimap2
        conf["minimap2"]["x_ava"] = "ava-pb"
        conf["minimap2"]["x_map"] = "asm20"
        # racon medaka
        conf["racon"]["iter"] = 2
        conf["medaka"]["iter"] = 0
        # UMI
        conf["umi"]["s"] = 60
        conf["umi"]["e"] = 60
        conf["umi"]["u"] = 3
        conf["umi"]["O"] = 0.05
        # simulate
        conf["simulate"]["badread"]["error_model"] = "pacbio2016"
        conf["simulate"]["badread"]["qscore_model"] = "pacbio2016"
        conf["simulate"]["badread"]["identity"] = "20,4"
        
    if longumi == True:
        conf["pool"] = False
        conf["primer_check"] = False
        conf["seqkit"]["min_qual"] = -1
        conf["seqkit"]["min_len"] = 3500
        conf["seqkit"]["max_len"] = 6000
        conf["fprimer"].clear()
        conf["fprimer"]["F"] = "AGRGTTYGATYMTGGCTCAG"
        conf["rprimer"].clear()
        conf["rprimer"]["R"] = "CGACATCGAGGTGCCAAAC"
        conf["fprimer_max"].clear()
        conf["fprimer_max"]["F"] = "AGRGTTYGATYMTGGCTCAG"
        conf["rprimer_min"].clear()
        conf["rprimer_min"]["R"] = "CGACATCGAGGTGCCAAAC"
        conf["umi"]["len"] = 18
        conf["umi"]["pattern"] = "NNNYRNNNYRNNNYRNNN NNNYRNNNYRNNNYRNNN"
        conf["umi"]["cutadapt"]["min_overlap"] = 11
        conf["flinker"] = "CAAGCAGAAGACGGCATACGAGAT"
        conf["rlinker"] = "AATGATACGGCGACCACCGAGATC"
        conf["min_cluster_size"] = 3
        conf["min_support_reads"] = 3
        
    if simulate == True:
        conf["fprimer"].clear()
        conf["fprimer"]["F"] = "AATGTACTTCGTTCAGTTACGTATTGCT"
        conf["rprimer"].clear()
        conf["rprimer"]["R"] = "ACTTCGTTCAGTTACGTATTGC"
        conf["fprimer_max"].clear()
        conf["fprimer_max"]["F"] = "AATGTACTTCGTTCAGTTACGTATTGCT"
        conf["rprimer_min"].clear()
        conf["rprimer_min"]["R"] = "ACTTCGTTCAGTTACGTATTGC"
        # simulate
        conf["simulate"]["badread"]["start_adapter_seq"] = "AATGTACTTCGTTCAGTTACGTATTGCT"
        conf["simulate"]["badread"]["end_adapter_seq"] = "GCAATACGTAACTGAACGAAGT"
        # medaka 
        conf["medaka"]["iter"] = 0
        # cluster
        conf["umapclust"]["hdbscan"]["min_bin_size"] = 10
        conf["umapclust"]["hdbscan"]["min_samples"] = 10
        conf["min_cluster_size"] = 10
            
    conf["basecalled_dir"] = bascdir
    conf["demultiplexed_dir"] = demuxdir
    # combine LACA runs, append if not empty
    merge_runs = []
    if merge_parent is not None:
        for parent in merge_parent:
            for sub in os.listdir(parent):
                if os.path.isdir(os.path.join(parent, sub)):
                    merge_runs.append(os.path.join(parent, sub)) 
        merge_runs += list(merge)
    else:
        merge_runs = list(merge)
    # remove duplicates and sort
    if len(merge_runs) > 0:
        conf["merge_runs"] = sorted(list(set(merge_runs)))
        
    conf["database_dir"] = dbdir
    conf["demuxer"] = demuxer
    conf["nreads_m"] = nreads_m
    conf["subsample"] = subsample
    conf["cluster"] = list(cluster)
    conf["consensus"] = list(consensus)
    conf["quant"] = list(quant)
    conf["uchime"] = uchime
    conf["threads"]["normal"] = jobs_m
    conf["threads"]["large"] = jobs_M
    
    if os.path.exists(os.path.join(workdir, config)):
        logger.warning(f"Config file [{config}] already exists in {workdir}.")
    else:
        with open(os.path.join(workdir, config), "w") as conf_file:
            yaml.dump(conf, conf_file)
        logger.info(f"Config file [{config}] created in {workdir}.")
    