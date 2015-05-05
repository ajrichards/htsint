************
Introduction
************

:Version: 0.4.3
:Authors: Adam Richards
:Web site: https://github.com/ajrichards/hts-integrate
:Documentation: http://ajrichards.github.io/hts-integrate
:Copyright: This document has been placed in the public domain.
:License: HTS-Integrate is released under the MIT License


Purpose
=======

``hts-integrate`` is a Python package principally used to create gene sets that are used in RNA-Seq studies.  Annotation data are integrated through the creation of gene sets.  The main functions included in ``hts-integrate`` are: 

  * Name mapping
  * BLAST
  * Gene Ontology queries
  * Heatmaps for differential expression analysis
  * Creation of gene sets for gene set enrichment analysis
  * Visualization of gene sets


Installation
================

For more details visit the documentation:

  *  http://ajrichards.github.io/hts-integrate

The easiest way to install and maintain ``hts-integrate`` is to use `pip <https://pypi.python.org/pypi/pip>`_

  .. code-block:: bash

      ~$ pip install htsint

  .. code-block:: bash

      ~$ pip install htsint --upgrade

Next you need to `setup your database <http://ajrichards.github.io/hts-integrate/database.html>`_ to get a working version of ``hts-integrate``.

Prerequsites 
-----------------------------------

  * `PostgreSQL <www.postgresql.org/>`_
  * https://networkx.github.io/
  * http://www.sqlalchemy.org

