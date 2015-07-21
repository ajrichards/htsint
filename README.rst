************
Introduction
************

:Version: 0.5.1
:Authors: Adam Richards
:Web site: https://github.com/ajrichards/htsint
:Documentation: http://ajrichards.github.io/htsint
:Copyright: This document has been placed in the public domain.
:License: htsint is released under the MIT License


Purpose
=======

``htsint`` (High-Throughput Sequencing INTegrate) is a Python package used to create gene sets for the study of high-throughput sequencing data. The goal is to create functional modules through the integration of heterogeneous types of data. These functional modules are primarily based on the Gene Ontology, but as the package matures additional sources of data will be incorporated. The functional modules produced can be subsequently tested for significance in terms of differential expression in RNA-Seq or microarray studies using gene set enrichment analysis. Shown below is the placement of htsint in an example analysis pipeline.

  * BLAST mapping
  * Gene Ontology queries
  * Heatmaps for differential expression analysis
  * Creation of gene sets for gene set enrichment analysis
  * Visualization of gene sets

Installation
================

For more details visit the documentation:

  *  http://ajrichards.github.io/htsint

The easiest way to install and maintain ``htsint`` is to use `pip <https://pypi.python.org/pypi/pip>`_

  .. code-block:: bash

      ~$ pip install htsint

  .. code-block:: bash

      ~$ pip install htsint --upgrade

Next you need to `setup your database <http://ajrichards.github.io/htsint/database.html>`_ to get a working version of ``htsint``.

Prerequsites 
-----------------------------------

  * `PostgreSQL <www.postgresql.org/>`_
  * `NetworkX <https://networkx.github.io/>`_
  * `SQLAlchemy <http://www.sqlalchemy.org/>`_
  * `Psycopg2 <http://initd.org/psycopg/>`_
  * `NumPy <www.numpy.org/>`_

  
