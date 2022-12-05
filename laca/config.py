from .log import logger

import os
from ruamel.yaml import YAML

def init_conf(
    bascdir,
    demuxdir,
    merge,
    dbdir,
    workdir,
    config="config.yaml",
    demuxer = "guppy",
    nreads_m = 1000,
    no_pool=False,
    subsample=False,
    no_trim=False,
    kmerbin=False,
    cluster=["isONclustCon"],
    chimer_filter=False,
    jobs_m=2,
    jobs_M=6,
    nanopore=False,
    pacbio=False,
    longumi=False,
):
    """
    Reads template config file with comments from ./template_config.yaml
    updates it by the parameters provided.
    Args:
        bascdir (str): path to a directory of basecalled fastq files
        demuxdir (str): path to a directory of demultiplexed fastq files
        merge (list): list of the working directories of LACA runs to merge
        dbdir (str): path to the taxonomy database
        workdir (str): path to the working directory
        config (str): the config filename
        demuxer (str): the demultiplexer [default: "guppy"]
        nreads_m (int): minimum number of reads for the demultiplexed fastqs
        no_pool (bool): if True, do not pool the reads [default: False]
        subsample (bool): if True, subsample the reads [default: False]
        no_trim (bool): if True, do not trim the primers [default: False]
        kmerbin (bool): if True, conduct kmer binning  [default: False]
        cluster (list): list of methods to generate consensus (kmerCon, clustCon, isONclustCon, isONcorCon, umiCon) [default: "isONclustCon"]
        chimer_filter (bool): if True, filter possible chimeras by vsearch [default: False]
        jobs_m (int): number of jobs for common tasks [default: 2]
        jobs_M (int): number of jobs for threads-dependent tasks [default: 6]
        nanopore (bool): if True, use template for nanopore reads [default: False]
        pacbio (bool): if True, use template for pacbio reads [default: False]
        longumi (bool): if True, use primer design from longumi paper (https://doi.org/10.1038/s41592-020-01041-y) [default: False]       
    """
    os.makedirs(dbdir, exist_ok=True)
    os.makedirs(workdir, exist_ok=True)
    
    yaml = YAML()
    template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workflow/config_template.yaml")
    
    with open(template_conf_file) as template_conf:
        conf = yaml.load(template_conf)
        
    if nanopore == True:
        # isONclust
        conf["isONclust"]["k"] = 13
        conf["isONclust"]["w"] = 20
        # minimap2
        conf["minimap2"]["x_ava"] = "ava-ont"
        conf["minimap2"]["x_map"] = "map-ont"
        # racon medaka
        conf["racon"]["iter"] = 2
        conf["medaka"]["iter"] = 1
        # isONcor
        conf["isONcor"] = True
        # UMI
        conf["umi"]["s"] = 90
        conf["umi"]["e"] = 90
        conf["umi"]["u"] = 3.5
        conf["umi"]["O"] = 0.2
    
    if pacbio == True:
        # isONclust
        conf["isONclust"]["k"] = 15
        conf["isONclust"]["w"] = 50
        # minimap2
        conf["minimap2"]["x_ava"] = "ava-pb"
        conf["minimap2"]["x_map"] = "asm20"
        # racon medaka
        conf["racon"]["iter"] = 2
        conf["medaka"]["iter"] = 0
        # isONcor
        conf["isONcor"] = False
        # UMI
        conf["umi"]["s"] = 60
        conf["umi"]["e"] = 60
        conf["umi"]["u"] = 3
        conf["umi"]["O"] = 0.05
        
    if longumi == True:
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
        conf["umi"]["seqkit"]["min_qual"] = -1
        conf["umi"]["seqkit"]["min_len"] = 3500
        conf["umi"]["seqkit"]["max_len"] = 6000
        conf["umi"]["len"] = 18
        conf["umi"]["pattern"] = "NNNYRNNNYRNNNYRNNN NNNYRNNNYRNNNYRNNN"
        conf["umi"]["cutadapt"]["min_overlap"] = 11
        conf["flinker"] = "CAAGCAGAAGACGGCATACGAGAT"
        conf["rlinker"] = "AATGATACGGCGACCACCGAGATC"
            
    conf["basecalled_dir"] = bascdir
    conf["demultiplexed_dir"] = demuxdir
    conf["merge_runs"] = list(merge)
    conf["database_dir"] = dbdir
    conf["demuxer"] = demuxer
    conf["nreads_m"] = nreads_m
    conf["pool"] = not no_pool
    conf["subsample"] = subsample
    conf["trim"] = not no_trim
    conf["kmerbin"] = kmerbin
    conf["cluster"] = list(cluster)
    conf["chimeraF"] = chimer_filter
    conf["threads"]["normal"] = jobs_m
    conf["threads"]["large"] = jobs_M
    
    if os.path.exists(os.path.join(workdir, config)):
        logger.warning(f"Config file [{config}] already exists in {workdir}.")
    else:
        with open(os.path.join(workdir, config), "w") as conf_file:
            yaml.dump(conf, conf_file)
        logger.info(f"Config file [{config}] created in {workdir}.")
    