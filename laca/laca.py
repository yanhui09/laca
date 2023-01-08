from email.policy import default
import click
import os
import subprocess
from .log import logger
from .config import init_conf

from snakemake import load_configfile
from . import __version__

def handle_max_mem(max_mem, profile):
    "Specify maximum memory (GB) to use. Memory is controlled by profile in cluster execution."
    "For numbers >1 its the memory in GB. "
    "For numbers <1 it's the fraction of available memory."

    if profile is not None:

        if max_mem is not None:
            logger.info(
                "Memory requirements are handled by the profile, I ignore max-mem argument."
            )
        # memory is handled via the profile, user should know what he is doing
        return ""
    else:
        import psutil
        from math import floor

        # calulate max  system meory in GB (float!)
        max_system_memory = psutil.virtual_memory().total / (1024**3)

        if max_mem is None:
            max_mem = 0.95
        if max_mem > 1:

            if max_mem > max_system_memory:
                logger.critical(
                    f"You specified {max_mem} GB as maximum memory, but your system only has {floor(max_system_memory)} GB"
                )
                sys.exit(1)

        else:

            max_mem = max_mem * max_system_memory

        # specify max_mem_string including java mem and max mem

        return f" --resources mem={floor(max_mem)} mem_mb={floor(max_mem*1024)} java_mem={floor(0.85* max_mem)} "

def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def run_smk(workflow, workdir, configfile, jobs, maxmem, profile, dryrun, snake_args, snakefile, exit_on_error, suppress):
    """
    Run Long Amplicon Consensus analysis (LACA).
    Most snakemake arguments can be appended, for more info see 'snakemake --help'
    """
    if not suppress:
        logger.info(f"LACA version: {__version__}")
    if not os.path.exists(configfile):
        logger.critical(f"Config file not found: {configfile}\nGenerate a config file using 'laca init'")
        exit(1)
    
    conf = load_configfile(configfile)
    # if snake_args (tuple, whitespace removed) contain "basecalled_dir=", take everything after as basecalled_dir
    for arg in snake_args:
        if "basecalled_dir=" in arg:
            conf["basecalled_dir"] = arg.split("=")[1]
    basecalled_dir = conf["basecalled_dir"]
    db_dir = conf["database_dir"]
    cmd = (
        "snakemake "
        "{wf} "
        "--directory '{workdir}' "
        "--snakefile '{snakefile}' "
        "--configfile '{configfile}' "
        "--use-conda {conda_prefix} "
        "{singularity_prefix} "
        "{singularity_args} "
        "{dryrun} "
        "--rerun-triggers mtime --rerun-incomplete --scheduler greedy "
        "--jobs {jobs} --nolock "
        " {max_mem} "
        " {profile} "
        " {args} "
    ).format(
        wf=workflow if workflow is not None else "",
        workdir=workdir,
        snakefile=snakefile,
        configfile=configfile,
        conda_prefix="--conda-prefix '" + os.path.join(db_dir, "conda_envs") + "'",
        singularity_prefix="--use-singularity --singularity-prefix '" + os.path.join(db_dir, "singularity_envs") + "'"
        if basecalled_dir is not None else "",
        singularity_args="--singularity-args '--bind " +
        os.path.dirname(snakefile) + "/resources/guppy_barcoding/:/opt/ont/guppy/data/barcoding/," + basecalled_dir + "'"
        if basecalled_dir is not None else "",
        dryrun="--dryrun" if dryrun else "",
        jobs=int(jobs) if jobs is not None else 1,
        max_mem=handle_max_mem(maxmem, profile),
        profile="" if (profile is None) else "--profile {}".format(profile),
        args=" ".join(snake_args),
    )
    if not suppress:
        logger.debug("Executing: %s" % cmd)
    
    try:
        if suppress:
            subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        if exit_on_error is True:
            exit(1)

# custom alo(at least one) class with mutex (mutually exclusive) classs
# https://stackoverflow.com/questions/44247099/click-command-line-interfaces-make-options-required-if-other-optional-option-is/
class AloMutex(click.Option):
    def __init__(self, *args, **kwargs):
        self.required_if_not:list = kwargs.pop("required_if_not")
        self.not_required_if:list = kwargs.pop("not_required_if")

        if self.required_if_not:
            kwargs["help"] = (
                kwargs.get("help", "") + " Option is required if '" + "', '".join(self.required_if_not) + "' not provided."
                ).strip()
        if self.not_required_if:
            kwargs["help"] = (
                kwargs.get("help", "") + " Option is mutually exclusive with '" + "', '".join(self.not_required_if) + "'."
                ).strip()
        super(AloMutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt:bool = self.name in opts
        if self.required_if_not and all(alo_opt not in opts for alo_opt in self.required_if_not) and not current_opt:
            raise click.UsageError(
                "at least one of '" + "', '".join(self.required_if_not) + 
                "' and '" + str(self.name) +  "' options is provided."
                )
        if self.not_required_if:
            for mutex_opt in self.not_required_if:
                if mutex_opt in opts:
                    if current_opt:
                        raise click.UsageError(
                            "'" + str(self.name) + "' is mutually exclusive with '" + str(mutex_opt) + "'."
                            )
                    else:
                        self.prompt = None
        return super(AloMutex, self).handle_parse_result(ctx, opts, args)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(
    __version__,
    "-v",
    "--version",
    )
@click.pass_context
def cli(self):
    """
    LACA: a reproducible and scaleable workflow for Long Amplicon Consensus Analysis.
    To follow updates and report issues, see: https://github.com/yanhui09/laca.
    """
    pass

@cli.command(
    'run',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run LACA workflow.'
)
@click.option(
    "-w",
    "--workdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Working directory for LACA.",
    show_default=True,
    default=".",
)
@click.option(
    "-c",
    "--configfile",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Config file for LACA. Use config.yaml in working directory if not specified.",
    default=None,
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=6,
    show_default=True,
    help="Maximum jobs to run in parallel.",
)
@click.option(
    "-m",
    "--maxmem",
    type=float,
    default=None,
    show_default=True,
    help=handle_max_mem.__doc__,
)
@click.option(
    "--profile",
    default=None,
    help="Snakemake profile for cluster execution.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Dry run.",
)

@click.argument(
    "workflow",
    type=click.Choice(
        ["demux", "qc", "kmerBin",
         "kmerCon", "clustCon", "isONclustCon", "isONclustCon2","isONcorCon", "umiCon",
         "quant", "taxa", "tree", "all", "merge",
         "initDB", "simulate"]
    ),
)
@click.argument("snake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(workflow, workdir, configfile, jobs, maxmem, profile, dryrun, snake_args):
    """
    Run LACA workflow.
    """
    if workflow == "simulate":
        sf = "workflow/rules/simulate.smk"
    else:
        sf = "workflow/Snakefile" 
    snakefile = get_snakefile(sf)
    configfile_run = os.path.join(workdir, "config.yaml") if configfile is None else configfile
    run_smk(workflow, workdir, configfile_run, jobs, maxmem, profile, dryrun, snake_args, snakefile, exit_on_error=True, suppress=False)

# laca init
# initialize config file
@cli.command(
    'init',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Prepare the config file.',
)
@click.option(
    "-b",
    "--bascdir",
    help="Path to a directory of the basecalled fastq files.",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["merge", "merge_parent", "demuxdir"],
)
@click.option(
    "-x",
    "--demuxdir",
    help="Path to a directory of demultiplexed fastq files.",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["merge", "merge_parent", "bascdir"],
)
@click.option(
    "--merge",
    help="Path to the working directory of a completed LACA run  [Mutiple]. Runs will be combined if --merge_parent applied.",
    multiple=True,
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["bascdir", "demuxdir"],
)
@click.option(
    "--merge-parent",
    help="Path to the parent of the working directories of completed LACA runs. Runs will be combined if --merge applied.",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["bascdir", "demuxdir"],
)
@click.option(
    "-d",
    "--dbdir",
    help="Path to the taxonomy databases.",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    "-w",
    "--workdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Output directory for LACA.",
    show_default=True,
    default=".",
)
@click.option(
    "--demuxer",
    type=click.Choice(["guppy", "minibar"]),
    default="guppy",
    show_default=True,
    help="Demultiplexer.",
)
@click.option(
    "--fqs-min",
    type=int,
    default=1000,
    show_default=True,
    help="Minimum number of reads for the demultiplexed fastqs.",
)
@click.option(
    "--no-pool",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not pool the reads for denoising.",
)
@click.option(
    "--subsample",
    is_flag=True,
    default=False,
    show_default=True,
    help="Subsample the reads.",
)
@click.option(
    "--no-trim",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not trim the primers.",
)
@click.option(
    "--kmerbin",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do pre-read binning.",
)
@click.option(
    "--cluster",
    type=click.Choice(
        ["kmerCon", "clustCon", "isONclustCon", "isONclustCon2", "isONcorCon", "umiCon"]
        ),
    default=["isONclustCon"],
    show_default=True,
    multiple=True,
    help="Consensus methods.  [Mutiple]",
)
@click.option(
    "--quant",
    type=click.Choice(["seqid", "minimap2"]),
    default=["seqid"],
    show_default=True,
    multiple=True,
    help="Create abundance matrix by sequence id or minimap2.  [Mutiple]",
)
@click.option(
    "--uchime",
    is_flag=True,
    default=False,
    show_default=True,
    help="Filter chimeras by uchime-denovo in vsearch.",
)
@click.option(
    "--jobs-min",
    type=int,
    default=2,
    show_default=True,
    help="Number of jobs for common tasks.",
)
@click.option(
    "--jobs-max",
    type=int,
    default=6,
    show_default=True,
    help="Number of jobs for threads-dependent tasks.",
)
@click.option(
    "--nanopore",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use config template for nanopore reads.",
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["isoseq"],
)
@click.option(
    "--isoseq",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use config template for pacbio CCS reads.",
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["nanopore"],
)
@click.option(
    "--longumi",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use primer design from longumi paper (https://doi.org/10.1038/s41592-020-01041-y).",
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["simulate"],
)
@click.option(
    "--simulate",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use config template for in silicon test.",
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["longumi"],
)
@click.option(
    "--clean-flags",
    is_flag=True,
    default=False,
    show_default=True,
    help="Clean flag files.",
)
def run_init(
    bascdir, demuxdir, merge, merge_parent, dbdir, workdir, demuxer, fqs_min, no_pool, subsample, no_trim, 
    kmerbin, cluster, quant, uchime, jobs_min, jobs_max, nanopore, isoseq, longumi, simulate, clean_flags):
    """
    Prepare config file for LACA.
    """ 
    logger.info(f"LACA version: {__version__}")
    init_conf(
        bascdir, demuxdir, merge, merge_parent, dbdir, workdir, "config.yaml", demuxer, fqs_min, no_pool, subsample,
        no_trim, kmerbin, cluster, quant, uchime, jobs_min, jobs_max, nanopore, isoseq, longumi, simulate)
    # clean flags if requested
    if clean_flags:
        # rm .*_DONE in workdir
        flags = [".qc_DONE"]
        for flag in flags:
            file = os.path.join(workdir, flag)
            if os.path.exists(file):
                os.remove(file)
        logger.warning(f"All flags files (.*_DONE) in {workdir} are removed.")
    
if __name__ == "__main__":
    cli()
