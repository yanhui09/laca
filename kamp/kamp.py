import click
import os
import logging
import subprocess

from snakemake import load_configfile
from snakemake.utils import validate
from .__init__ import __version__


def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

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
    short_help='run kamp workflow'
)
@click.option(
    "-c",
    "--configfile",
    type=click.Path(exists=True, resolve_path=True),
    help="Config file to run kamp.",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=1,
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
         "nanosim"]
    ),
)
@click.argument("snake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(
    workflow, configfile, jobs, maxmem, dryrun, snake_args
):
    """
    Run Kamp workflow to proceess long read amplicons.
    Most snakemake arguments can be appended, for more info see 'snakemake --help'
    """
    
    logger = logging.getLogger()
    logger.info(f"Kamp version: {__version__}")

    if not os.path.exists(configfile):
        logger.critical(
            f"config-file not found: {configfile}\n"
        )
        exit(1)
    
    conf = load_configfile(configfile)

    db_dir = conf["database_dir"]

    cmd = (
        "snakemake {wf} --snakefile {snakefile} "
        "--configfile '{configfile}' "
        "--use-conda {conda_prefix} {dryrun} "
        "--rerun-triggers mtime --rerun-incomplete "
        "--jobs {jobs} --nolock "
        " {max_mem} "
        " {args} "
    ).format(
        wf=workflow if workflow != "None" else "",
        snakefile=get_snakefile(),
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

if __name__ == "__main__":
    cli()