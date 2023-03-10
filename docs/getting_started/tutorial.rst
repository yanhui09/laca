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

.. code-block:: bash

    laca run all -j 50 -m 100

Cluster execuation
------------------

.. code-block:: bash

    laca run all --profile cluster -j 20

Use ``snakemake`` arguments
---------------------------

``LACA`` wraps ``snakemake`` to excuate the workflow, and accepts ``--snakemake`` arguments 
in execuation. Except ``--profile``, other useful arguments includes:

- ``-n``, which prints the dry-run overview.
- ``-k/--keep-going``, which allows laca run to acomplish other independent jobs even if some jobs fail.
- ``--report``, which generates a report of the workflow by specifying e.g. ``--report report.html``.

More details can be found in the `snakemake documentation <https://snakemake.readthedocs.io/en/stable/>`_.

Merge runs
==========

``LACA`` allows meta-analysis of multiple accomplished runs, and generates the integrated OTU catalogue, count matrix 
and annotations. To merge runs, you can use the flag ``-m/--merge`` in ::ref::`laca init` to specify the path to the 
accomplished runs.

.. code-block:: bash

    laca init -m /path/to/laca_run1 -m /path/to/laca_run2 -d /path/to/database
    laca run merge

Run with simulated reads
========================

``LACA`` supports performance benchmark between genetic heterogeneity of target region, sequencing noise and depth. 
The *in silico* long amplicon reads are generated with `Badread <https://github.com/rrwick/Badread>`_ 
from `SILVA SSUs <https://www.arb-silva.de/>`_ according the defined parameters in the ``config`` file. To enable this, 

.. code-block:: bash

    laca init -d /path/to/database --simulate
    laca run simulate
    laca run all
