.. primertool documentation master file, created by
   sphinx-quickstart on Tue Jul  9 12:33:04 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Primertool
======================================

The premise of this package is to generate primers for PCR/Sanger sequencing for either:

- a specific variant (hgvs_ variant nomenclature), either the whole exon or if not in an exon for the genomic position
- an exon (transcript number and exon )
- all exons of a transcript (transcript number)
- around a genomic position (chromosome and start/stop position)

This tool allows for primers based the human reference genome version hg19 or hg38.

This project can be used in three ways. The core functionality of this tool is contained in the primertool python
package inside this project and can be used via the command line (usage examples below). The streamlit app also
contained in this project is simply a wrapper around the core functionality of the primertool package to make the tool
accessible to non-programmers. The streamlit app can be run locally or deployed to a server. For proper deployment and
management a Dockerfile is also provided.

Usage options:

- python package
- streamlit app (recommended for development)
- docker container (recommended for deployment)

Overview
--------

For a detailed overview of how the tool generates primers, please refer to the :doc:`information <information>` page.
The software documentation is available in the :doc:`modules <modules>` page. For instructions on how to set up the
tool, please refer to the :doc:`setup <setup>` page.

Installation
--------------

Using the python package requires Python 3.8+. All required packages are listed in the requirements file.

Local installation is possible by running

.. code:: bash

    conda env create -f environment.yml

in the main directory (requires conda and pip installed).

You will also need gcc installed on your system, so check if you already have it by running

.. code:: bash

    gcc --version

If you don't have gcc installed, you can install it by running

.. code:: bash

    sudo apt update
    sudo apt install build-essential

For the streamlit app and docker container, the setup instructions are provided in the :doc:`setup guide <setup>`.
When running the tool behind a firewall and it is not possible to access the hg38 genome database (or a local version
is wanted for other reasons), it is possible to run a local instance of said database. To create such a local MySQL
instance follow the steps described here: :doc:`MySQL setup guide <mysql-setup>`

Basic usage examples
---------------------

1. Generating primers for a mutation using VariantPrimerGenerator:

.. code:: Python

    from primertool import primertool as pt
    df = pt.VariantPrimerGenerator('NM_003165.6:c.1702G>A', 'hg38', kuerzel='SM').ordertable

2. Generating primer for a single exon using ExonPrimerGenerator:

.. code:: Python

    from primertool import primertool as pt
    df = pt.ExonPrimerGenerator('NM_000451', 6, 'hg38', kuerzel='SM').ordertable

3. Generating primers for every exon in a gene using GenePrimerGenerator:

.. code:: Python

    from primertool import primertool as pt
    df = pt.GenePrimerGenerator('NM_000451', 'hg38', kuerzel='SM').ordertable

4. Generating primer for a genomic position using GenomicPositionPrimerGenerator:

.. code:: Python

    from primertool import primertool as pt
    df = pt.GenomicPositionPrimerGenerator('chr12', 121814175, 121814175, 'hg38', kuerzel='SM').ordertable

.. toctree::
   information
   :maxdepth: 2
   :caption: Contents

.. toctree::
   modules
   :maxdepth: 2
   :caption: Package

.. toctree::
   setup
   mysql-setup
   deployment-firewall
   :maxdepth: 2
   :caption: Setup and Deployment


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. Links
.. _hgvs: https://varnomen.hgvs.org/
