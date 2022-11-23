from email.policy import default
import click
import os
import subprocess
from .log import logger
from .config import init_conf

from snakemake import load_configfile
from .__init__ import __version__


def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def run_smk(workflow, workdir, configfile, jobs, maxmem, dryrun, snake_args, snakefile, exit_on_error, suppress):
    """
    Run Kamp workflow to proceess long read amplicons.
    Most snakemake arguments can be appended, for more info see 'snakemake --help'
    """
    if not suppress:
        logger.info(f"Kamp version: {__version__}")
    if not os.path.exists(configfile):
        logger.critical(f"Config file not found: {configfile}\nGenerate a config file using 'kamp init'")
        exit(1)
    
    conf = load_configfile(configfile)
    db_dir = conf["database_dir"]
    basecalled_dir = conf["basecalled_dir"]
    cmd = (
        "snakemake "
        "{wf} "
        "--directory {workdir} "
        "--snakefile {snakefile} "
        "--configfile '{configfile}' "
        "--use-conda {conda_prefix} "
        "{singularity_prefix} "
        "{singularity_args} "
        "{dryrun} "
        "--rerun-triggers mtime --rerun-incomplete --scheduler greedy "
        "--jobs {jobs} --nolock "
        " {max_mem} "
        " {args} "
    ).format(
        wf=workflow if workflow is not None else "",
        workdir=workdir,
        snakefile=snakefile,
        configfile=configfile,
        conda_prefix="--conda-prefix " + os.path.join(db_dir, "conda_envs"),
        singularity_prefix="--use-singularity --singularity-prefix " + os.path.join(db_dir, "singularity_envs")
        if basecalled_dir != None else "",
        singularity_args="--singularity-args '--bind " +
        os.path.dirname(snakefile) + "/resources/guppy_barcoding/:/opt/ont/guppy/data/barcoding/," + basecalled_dir + "'"
        if basecalled_dir is not None else "",
        dryrun="--dryrun" if dryrun else "",
        jobs=int(jobs) if jobs is not None else 1,
        max_mem="--resources mem_mb={}".format(int(float(maxmem)*1024)) if maxmem is not None else 50,
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

# custom alo(at least one) class from mutex (mutually exclusive) classs
# https://stackoverflow.com/questions/44247099/click-command-line-interfaces-make-options-required-if-other-optional-option-is/
class Alo(click.Option):
    def __init__(self, *args, **kwargs):
        self.required_if_not:list = kwargs.pop("required_if_not")

        assert self.required_if_not, "'required_if_not' parameter required"
        kwargs["help"] = (kwargs.get("help", "") + "Option is required if '" + ", ".join(self.required_if_not) + "' not provided.").strip()
        super(Alo, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt:bool = self.name in opts
        for alo_opt in self.required_if_not:
            if alo_opt not in opts:
                if not current_opt:
                    raise click.UsageError("At least one of '" + str(self.name) + "' or '" + str(alo_opt) + "' is required.")
                else:
                    self.prompt = None
        return super(Alo, self).handle_parse_result(ctx, opts, args)

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
    short_help='Run Kamp workflow.'
)
@click.option(
    "-w",
    "--workdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Output directory for Kamp.",
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
        ["demux", "qc", "kmerBin",
         "kmerCon", "clustCon", "isONclustCon", "isONcorCon", "umiCon",
         "quant", "taxa", "tree", "requant", "all",
         "initDB", "nanosim"]
    ),
)
@click.argument("snake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(workflow, workdir, jobs, maxmem, dryrun, snake_args):
    """
    Run Kamp workflow.
    """
    if workflow == "nanosim":
        sf = "workflow/rules/nanosim.smk"
    else:
        sf = "workflow/Snakefile" 
    snakefile = get_snakefile(sf)
    configfile = os.path.join(workdir, "config.yaml")
    run_smk(workflow, workdir, configfile, jobs, maxmem, dryrun, snake_args, snakefile, exit_on_error=True, suppress=False)

# kamp init
# initialize config file
@cli.command(
    'init',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Prepare the config file.',
)
@click.option(
    '-b',
    '--bascdir',
    help='Path to a directory of the basecalled fastq files. ',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = Alo,
    required_if_not = ['demuxdir'],
)
@click.option(
    '-x',
    '--demuxdir',
    help='Path to a directory of demultiplexed fastq files. ',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = Alo,
    required_if_not = ['bascdir'],
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
    help="Output directory for Kamp.",
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
    help="Use kmer binning.",
)
@click.option(
    "--cluster",
    type=click.Choice(
        ["kmerCon", "clustCon", "isONclustCon", "isONcorCon", "umiCon"]
        ),
    default=["isONclustCon"],
    show_default=True,
    multiple=True,
    help="Consensus methods.",
)
@click.option(
    "--chimerf",
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
    help="Config template for nanopore reads.",
)
@click.option(
    "--pacbio",
    is_flag=True,
    default=False,
    show_default=True,
    help="Config template for pacbio CCS reads.",
)
@click.option(
    "--longumi",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use primer design from longumi paper (https://doi.org/10.1038/s41592-020-01041-y).",
)
@click.option(
    "--clean-flags",
    is_flag=True,
    default=False,
    show_default=True,
    help="Clean flag files.",
)
def run_init(
    bascdir, demuxdir, dbdir, workdir, demuxer, fqs_min, no_pool, subsample, no_trim, 
    kmerbin, cluster, chimerf, jobs_min, jobs_max, nanopore, pacbio, longumi, clean_flags):
    """
    Prepare config file for Kamp.
    """ 
    logger.info(f"Kamp version: {__version__}")
    init_conf(
        bascdir, demuxdir, dbdir, workdir, "config.yaml", demuxer, fqs_min, no_pool, subsample,
        no_trim, kmerbin, cluster, chimerf, jobs_min, jobs_max, nanopore, pacbio, longumi)
    # clean flags if requested
    if clean_flags:
        # rm .*_DONE in workdir
        flags = [".simulated_DONE", ".qc_DONE"]
        for flag in flags:
            file = os.path.join(workdir, flag)
            if os.path.exists(file):
                os.remove(file)
        logger.warning(f"All flags files (.*_DONE) in {workdir} are removed.")
    
if __name__ == "__main__":
    cli()