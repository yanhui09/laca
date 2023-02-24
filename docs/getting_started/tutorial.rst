Tutorial
********

``LACA`` is easy to use. You can start a new analysis in two steps using ``laca init`` and ``laca run`` . 

.. code-block:: bash

    conda activate laca                                             # activate required environment 
    laca init -b /path/to/basecalled_fastqs -d /path/to/database    # init config file and check
    laca run all                                                    # start analysis

Initialization 
==============

:ref:`laca_init` will generate a config file in the working directory, 
which contains the necessary parameters to run ``LACA``. To initialize 
a new run, you need to provide at least following information:

-  ``-b/--basecalled``: path to the basecalled ``fastq`` or ``fastq.gz`` files

   or ``-x/--demuxdir``: path to the demultiplexed ``fastq`` files

   or ``--merge/--merge-parent``: path to the accomplished ``LACA`` runs

-  ``-d/--dbdir``: path to the database directory
-  ``-w/--workdir``: path to the working directory [default: current directory]

The generated config file was utilized for excuate the ``snakemake`` workflow. 
We suggest you review and adapt the config file, according to your primers, 
polishing procedures, etc. Below are some examples to initialize ``LACA``.

Input data
----------

**Basecalled data**

.. _guppy: https://community.nanoporetech.com/downloads/guppy/release_notes
.. _minibar: https://github.com/calacademy-research/minibar

For the basecalled Oxford Nanopore reads, ``LACA`` uses guppy_ as the default demultiplexer and also 
supports minibar_. To use ``minibar`` as the demutiplexer, you can append the flag ``--demuxer minibar``. 
To start demultiplexing with ``guppy``,

.. code-block:: bash

    laca init -b /path/to/basecalled_fastqs -d /path/to/database

**Demultiplexed data**

To start ``LACA`` with demultiplexed data, you can use ``-x/--demuxdir`` to specify the directory path. 
The demultiplex directory follows specific structure:

- Demultiplexed ``fastq`` files are placed in independent subdirectories, 
  each of which is named by the ``barcode`` name.
- Only subdirectory name follows the format ``[A-Za-z]+[0-9]+`` will be recognized as ``barcode``. 

.. code-block:: bash

    laca init -x /path/to/demultiplexed_fastqs -d /path/to/database

Predefined templates
--------------------

The predefined templates can be used with appended flag  ``--nanopore``, ``--isoseq``,
``--longumi`` and ``--simulate``, respectively. To use the ``isoseq`` template for 
Pacbio CCS reads, 

.. code-block:: bash

    laca init -b /path/to/basecalled_fastqs -d /path/to/database --isoseq

Execuation
==========

:ref:`laca_run` will trigger the full workflow or a specfic module under defined resource accordingly.

Command Line interface
----------------------

Cluster execuation
------------------

Useful ``run`` arguments
------------------------

Get a dry-run overview with ``-n``. ``snakemake`` arguments can be appened to ``laca run`` as well.

Merge runs
==========

Read simulation
===============