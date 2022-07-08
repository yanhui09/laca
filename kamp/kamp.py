from email.policy import default
import click
import os
import subprocess
from .log import logger
from .config import init_conf

from snakemake import load_configfile
#from snakemake.utils import validate
from .__init__ import __version__


def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def run_smk(
    workflow, workdir, configfile, jobs, maxmem, dryrun, snake_args, snakefile,
):
    """
    Run Kamp workflow to proceess long read amplicons.
    Most snakemake arguments can be appended, for more info see 'snakemake --help'
    """
    
    logger.info(f"Kamp version: {__version__}")
    
    if not os.path.exists(configfile):
        logger.critical(
            f"Config file not found: {configfile}\nGenerate a config file using 'kamp init'"
        )
        exit(1)
    
    conf = load_configfile(configfile)

    db_dir = conf["database_dir"]

    cmd = (
        "snakemake "
        "{wf} "
        "--directory {workdir} "
        "--snakefile {snakefile} "
        "--configfile '{configfile}' "
        "--use-conda {conda_prefix} "
        "{dryrun} "
        "--rerun-triggers mtime --rerun-incomplete "
        "--jobs {jobs} --nolock "
        " {max_mem} "
        " {args} "
    ).format(
        wf=workflow if workflow != "None" else "",
        workdir=workdir,
        snakefile=snakefile,
        configfile=configfile,
        conda_prefix="--conda-prefix " + os.path.join(db_dir, "conda_envs"),
        dryrun="--dryrun" if dryrun else "",
        jobs=int(jobs) if jobs is not None else 1,
        max_mem="--resources mem_mb={}".format(int(float(maxmem)*1024)) if maxmem is not None else 50,
        args=" ".join(snake_args),
    )
    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        exit(1)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(
    __version__,
    "-v",
    "--version",
    )
@click.pass_context
def cli(self):
    """
    Kamp: a k-mer based denoise pipeline to process long read amplicon sequencing by Nanopore.
    To follow updates and report issues, see: https://github.com/yanhui09/Kamp.
    """
    pass

@cli.command(
    'run',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run kamp workflow'
)
@click.option(
    "-w",
    "--workdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Output directory for kamp.",
    show_default=True,
    default=".",
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
    default=50,
    show_default=True,
    help="Maximum memory to use in GB.",
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
        ["demultiplex", "qc", "kmerBin",
         "kmerCon", "clustCon", "isONclustCon", "isONcorCon", "umiCon",
         "quant", "taxa", "tree", "requant", "all",
         "initDB", "nanosim"]
    ),
)
@click.argument("snake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(workflow, workdir, jobs, maxmem, dryrun, snake_args):
    """
    Run the kamp main workflow.
    """
    if workflow == "nanosim":
        sf = "workflow/rules/nanosim.smk"
    else:
        sf = "workflow/Snakefile" 
    snakefile = get_snakefile(sf)
    configfile = os.path.join(workdir, "config.yaml")
    run_smk(
        workflow, workdir, configfile, jobs, maxmem, dryrun, snake_args, snakefile
    )

# kamp init
# initialize taxonomy database and config file
@cli.command(
    'init',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Prepare config file and taxonomy database.',
)
@click.option(
    '-f',
    '--fqdir',
    help='Path to the basecalled fastq files.',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    '-d',
    '--dbdir',
    help='Path to the taxonomy databases.',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    "-w",
    "--workdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Output directory for kamp.",
    show_default=True,
    default=".",
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
    help="Conduct kmer binning.",
)
@click.option(
    "--cluster",
    type=click.Choice(
        ["kmerCon", "clustCon", "isONclustCon", "isONcorCon", "umiCon"]
        ),
    default=["isONclustCon"],
    show_default=True,
    multiple=True,
    help="Methods to generate consensus.",
)
@click.option(
    "--chimerf",
    is_flag=True,
    default=False,
    show_default=True,
    help="Filter possible chimeras by vsearch.",
)
@click.option(
    "--jobs-min",
    type=int,
    default=2,
    show_default=True,
    help="Number of jobs for common tasks",
)
@click.option(
    "--jobs-max",
    type=int,
    default=6,
    show_default=True,
    help="Number of jobs for threads-dependent tasks",
)
def run_init(fqdir, dbdir, workdir, no_pool, subsample, no_trim, kmerbin, cluster, chimerf, jobs_min, jobs_max):
    """
    Prepare the config file and working directory.
    """ 
    init_conf(fqdir, dbdir, workdir, "config.yaml", 
              no_pool, subsample, no_trim, kmerbin, cluster, chimerf, jobs_min, jobs_max)

if __name__ == "__main__":
    cli()