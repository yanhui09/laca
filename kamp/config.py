from .log import logger

import os
from ruamel.yaml import YAML

def init_conf(
    fqdir,
    dbdir,
    workdir,
    config="config.yaml",
    demultiplex = "guppy",
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
):
    """
    Reads template config file with comments from ./template_config.yaml
    updates it by the parameters provided.
    Args:
        fqdir (str): path to the basecalled fastq files
        dbdir (str): path to the taxonomy database
        workdir (str): path to the working directory
        config (str): the config filename
        demultiplex (str): the demultiplexing method [default: "guppy"]
        nreads_m (int): minimum number of reads for the demultiplexed fastqs
        no_pool (bool): if True, do not pool the reads [default: False]
        subsample (bool): if True, subsample the reads [default: False]
        no_trim (bool): if True, do not trim the primers [default: False]
        kmerbin (bool): if True, conduct kmer binning  [default: False]
        cluster (str): list of methods to generate consensus (kmerCon, clustCon, isONclustCon, isONcorCon, umiCon) [default: "isONclustCon"]
        chimer_filter (bool): if True, filter possible chimeras by vsearch [default: False]
        jobs_m (int): number of jobs for common tasks [default: 2]
        jobs_M (int): number of jobs for threads-dependent tasks [default: 6]
        nanopore (bool): if True, use template for nanopore reads [default: False]
        pacbio (bool): if True, use template for pacbio reads [default: False]       
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
            
    conf["basecalled_dir"] = fqdir
    conf["database_dir"] = dbdir
    conf["demultiplex"] = demultiplex
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
    