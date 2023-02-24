Expected output
***************

all
===
.. code-block:: bash

  laca run all

It runs the whole workflow and generates:

  - ``count_matrix.tsv`` 
  - ``rep_seqs.fasta``
  - ``taxonomy.tsv`` and ``taxonomy/{classifier}/taxonomy.tsv`` per defined classifier
  - ``tree.nwk`` and ``tree/{phylogeny}/tree.nwk`` per defined phylogeny tree maker.

demux
=====
.. code-block:: bash

  laca run demux

Per barcode it generates:

  - ``demultiplexed/{barcode}``
  - ``qc/{barcode}.fastq``

qc
===
.. code-block:: bash

  laca run qc

Per barcode it generates:

  - ``qc/qfilt/{barcode}.fastq``

kmerBin
=======
.. code-block:: bash

  laca run binning

Per barcode it generates:

  - ``kmerBin/split/{barcode}_{kmerbin}_0.fastq`` if ``kmerbin=true`` in ``config.yaml``
  - ``qc/qfilt/{barcode}.fastq`` if ``kmerbin=false`` in ``config.yaml``

consensus
=========
.. code-block:: bash

  laca run kmerCon
  # or
  laca run clustCon
  # or
  laca run isONclustCon
  # or
  laca run isONclustCon2
  # or
  laca run isONcorCon
  # or
  laca run umiCon
  

``LACA`` supports multiple consensus-calling approaches and 
per approach it generates:

  - ``{consensus}/{consensus}.fna`` 

quant
=====
.. code-block:: bash

  laca run quant

According to the defined ``quant`` approach in ``config.yaml``, it generates:

  - ``count_matrix.tsv`` 
  - ``rep_seqs.fasta``

taxa
====
.. code-block:: bash

  laca run taxa

According to the defined ``classifier`` in ``config.yaml``, it generates:

  - ``taxonomy.tsv`` and ``taxonomy/{classifier}/taxonomy.tsv`` per defined classifier

tree
====
.. code-block:: bash

  laca run tree

According to the defined ``phylogeny`` tree maker in ``config.yaml``, it generates:

  - ``tree.nwk`` and ``tree/{phylogeny}/tree.nwk`` per defined phylogeny tree maker.

merge
=====
.. code-block:: bash


  laca run merge

According to the defined ``basecalled_dir``, ``demultiplexed_dir``, 
``merge_runs`` in ``config.yaml``, it generates:

  - ``count_matrix_merged.tsv`` and ``rep_seqs_merged.fasta``

  If ``basecalled_dir`` and ``demultiplexed_dir`` are ``None``, 
  and ``merge_runs`` is not ``None``, it re-generates the 
  taxonomy and tree files.

  - ``taxonomy.tsv`` and ``taxonomy/{classifier}/taxonomy.tsv`` per defined classifier
  - ``tree.nwk`` and ``tree/{phylogeny}/tree.nwk`` per defined phylogeny tree maker.
  
simulate
========
.. code-block:: bash

  laca run simulate

According to the ``config.yaml``, it prepares *in silico* demultiplexed reads:

  - ``demultiplexed/{barcode}``