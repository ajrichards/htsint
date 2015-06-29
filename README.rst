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

HTSInt (High-Througput Sequencing INTegrate) is a Python package principally used to create gene sets that are used in sequencing studies.  Annotation data are integrated through the creation of gene sets.  HTSInt also serves as a general purpose library to study gene expression.  The target audience for HTSInt are developers that piece together high-throughput sequencing (HTS) pipelines.  The software in its current form is an API library and because HTS pipelines have different goals with many varied tasks required to achieve these goals, a flexible library in a scripting language commonly used in bioinformatics was selected.  One important aspect of this project is to provide both abstracted functions for non-Python programmers as well as convenient means to enable a higher levels of customization.  The main functions included in HTSInt are: 

  * BLAST mapping
  * Gene Ontology queries
  * Heatmaps for differential expression analysis
  * Creation of gene sets for gene set enrichment analysis
  * Visualization of gene sets

Installation
================

For more details visit the documentation:

  *  http://ajrichards.github.io/htsint

The easiest way to install and maintain HTSInt is to use `pip <https://pypi.python.org/pypi/pip>`_

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

  
