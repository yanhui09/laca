User manual
***************

LACA Commands
=============

``LACA`` is easy to use. You can start a new analysis in two steps using ``laca init`` and ``laca run`` . 

.. code-block:: text

    Usage: laca [OPTIONS] COMMAND [ARGS]...
    
      LACA: a reproducible and scaleable workflow for Long Amplicon Consensus
      Analysis. To follow updates and report issues, see:
      https://github.com/yanhui09/laca.
    
    Options:
      -v, --version  Show the version and exit.
      -h, --help     Show this message and exit.
    
    Commands:
      init  Prepare the config file.
      run   Run LACA workflow.

LACA Init
=========

``laca init`` will generate a config file in the working directory, 
which contains the necessary parameters to run ``LACA``.

.. code-block:: text

    Usage: laca init [OPTIONS]
    
      Prepare config file for LACA.
    
    Options:
      -b, --bascdir PATH              Path to a directory of the basecalled fastq
                                      files. Option is mutually exclusive with
                                      'merge', 'merge_parent', 'demuxdir'.
      -x, --demuxdir PATH             Path to a directory of demultiplexed fastq
                                      files. Option is mutually exclusive with
                                      'merge', 'merge_parent', 'bascdir'.
      --merge PATH                    Path to the working directory of a completed
                                      LACA run  [Mutiple]. Runs will be combined
                                      if --merge_parent applied. Option is
                                      mutually exclusive with 'bascdir',
                                      'demuxdir'.
      --merge-parent PATH             Path to the parent of the working
                                      directories of completed LACA runs. Runs
                                      will be combined if --merge applied. Option
                                      is mutually exclusive with 'bascdir',
                                      'demuxdir'.
      -d, --dbdir PATH                Path to the taxonomy databases.  [required]
      -w, --workdir PATH              Output directory for LACA.  [default: .]
      --demuxer [guppy|minibar]       Demultiplexer.  [default: guppy]
      --fqs-min INTEGER               Minimum number of reads for the
                                      demultiplexed fastqs.  [default: 1000]
      --no-pool                       Do not pool the reads for denoising.
      --subsample                     Subsample the reads.
      --no-trim                       Do not trim the primers.
      --kmerbin                       Do pre-read binning.
      --cluster [kmerCon|clustCon|isONclustCon|isONclustCon2|isONcorCon|umiCon]
                                      Consensus methods.  [Mutiple]  [default:
                                      isONclustCon]
      --quant [seqid|minimap2]        Create abundance matrix by sequence id or
                                      minimap2.  [Mutiple]  [default: seqid]
      --uchime                        Filter chimeras by uchime-denovo in vsearch.
      --jobs-min INTEGER              Number of jobs for common tasks.  [default:
                                      2]
      --jobs-max INTEGER              Number of jobs for threads-dependent tasks.
                                      [default: 6]
      --nanopore                      Use config template for nanopore reads.
                                      Option is mutually exclusive with 'isoseq'.
      --isoseq                        Use config template for pacbio CCS reads.
                                      Option is mutually exclusive with
                                      'nanopore'.
      --longumi                       Use primer design from longumi paper (https:
                                      //doi.org/10.1038/s41592-020-01041-y).
                                      Option is mutually exclusive with
                                      'simulate'.
      --simulate                      Use config template for in silicon test.
                                      Option is mutually exclusive with 'longumi'.
      --clean-flags                   Clean flag files.
      -h, --help                      Show this message and exit.

LACA Run
========

``laca run`` will trigger the full workflow or a specfic module under defined resource accordingly.
Get a dry-run overview with ``-n``. ``snakemake`` arguments can be appened to ``laca run`` as well.

.. code-block:: text

    Usage: laca run [OPTIONS] {demux|qc|kmerBin|kmerCon|clustCon|isONclustCon|isON
                    clustCon2|isONcorCon|umiCon|quant|taxa|tree|all|merge|initDB|s
                    imulate} [SNAKE_ARGS]...
    
      Run LACA workflow.
    
    Options:
      -w, --workdir PATH     Working directory for LACA.  [default: .]
      -c, --configfile FILE  Config file for LACA. Use config.yaml in working
                             directory if not specified.
      -j, --jobs INTEGER     Maximum jobs to run in parallel.  [default: 6]
      -m, --maxmem FLOAT     Specify maximum memory (GB) to use. Memory is
                             controlled by profile in cluster execution.
      --profile TEXT         Snakemake profile for cluster execution.
      -n, --dryrun           Dry run.
      -h, --help             Show this message and exit.
