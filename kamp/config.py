from .log import logger

import os
from ruamel.yaml import YAML

def init_conf(
    fqdir,
    dbdir,
    workdir,
    config="config.yaml",
    no_pool=False,
    subsample=False,
    no_trim=False,
    kmerbin=False,
    cluster=["isONclustCon"],
    chimer_filter=False,
    jobs_m=2,
    jobs_M=6,
):
    """
    Reads template config file with comments from ./template_config.yaml
    updates it by the parameters provided.
    Args:
        fqdir (str): path to the basecalled fastq files
        dbdir (str): path to the taxonomy database
        workdir (str): path to the working directory
        config (str): the config filename
        no_pool (bool): if True, do not pool the reads [default: False]
        subsample (bool): if True, subsample the reads [default: False]
        no_trim (bool): if True, do not trim the primers [default: False]
        kmerbin (bool): if True, conduct kmer binning  [default: False]
        cluster (str): list of methods to generate consensus (kmerCon, clustCon, isONclustCon, isONcorCon, umiCon) [default: "isONclustCon"]
        chimer_filter (bool): if True, filter possible chimeras by vsearch [default: False]
        jobs_m (int): number of jobs for common tasks [default: 2]
        jobs_M (int): number of jobs for threads-dependent tasks [default: 6]       
    """
    os.makedirs(dbdir, exist_ok=True)
    os.makedirs(workdir, exist_ok=True)
    
    yaml = YAML()
    template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workflow/config.yaml")
    
    with open(template_conf_file) as template_conf:
        conf = yaml.load(template_conf)
        
    conf["basecalled_dir"] = fqdir
    conf["database_dir"] = dbdir
    conf["pool"] = not no_pool
    conf["subsample"] = subsample
    conf["trim_primers"] = not no_trim
    conf["kmerbin"] = kmerbin
    conf["cluster"] = list(cluster)
    conf["chimer_filter"] = chimer_filter
    conf["jobs_m"] = jobs_m
    conf["jobs_M"] = jobs_M
    
    if os.path.exists(os.path.join(workdir, config)):
        logger.warning(f"Config file [{config}] already exists in {workdir}.")
    else:
        with open(os.path.join(workdir, config), "w") as conf_file:
            yaml.dump(conf, conf_file)
        logger.info(f"Config file [{config}] created in {workdir}.")
    