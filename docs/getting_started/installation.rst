Installation
***************

Conda installation
==================

`Conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ is the only required dependency prior to installation.
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ is enough for whole installation. 

Install from GitHub
===================

1. Clone GitHub repository and create an independent ``conda`` environment

.. code-block:: bash

    git clone https://github.com/yanhui09/laca.git
    cd laca
    conda env create -n laca -f env.yaml 

Or speed it up with ``mamba``.

.. code-block:: bash

    mamba env create -n laca -f env.yaml 

2. Install ``LACA`` with ``pip``
      
To avoid inconsistency, we suggest installing ``LACA`` in the ``conda`` environment created previously.

.. code-block:: bash

    conda activate laca
    pip install --editable .

Installation check with demo
============================

.. code-block:: bash

    cd $HOME
    conda activate laca 
    laca init -b laca/workflow/resources/data -d test/Database -w test
    laca run all -w test -n 
