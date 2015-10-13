.. main file for htsint documentation

htsint
=====================================================

.. figure:: ../figures/example-web.png
   :scale: 55%
   :align: center
   :alt: example functional module
   :figclass: align-center

-------


``htsint`` (High-Throughput Sequencing INTegrate) is a Python package used to create gene sets for the study of high-throughput sequencing data.  The goal is to create functional modules through the integration of heterogeneous types of data.  These functional modules are primarily based on the Gene Ontology [Ashburner00]_, but as the package matures additional sources of data will be incorporated.  The functional modules produced can be subsequently tested for significance in terms of differential expression in RNA-Seq or microarray studies using gene set enrichment analysis [Subramian05]_.  Shown below is the placement of ``htsint`` in an example analysis pipeline.

.. figure:: ../figures/figure-conceptual.png
   :scale: 65%
   :align: center
   :alt: example RNA-Seq pipeline
   :figclass: align-center


Who are the intended users?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The software in its current form is an API library and because HTS pipelines have different goals with many varied tasks required to achieve these goals, a flexible library in a scripting language commonly used in bioinformatics was selected.  The target audience for ``htsint`` are developers that piece together high-throughput sequencing (HTS) pipelines.  That being said one important aspect of this project is to provide both abstracted functions for non-Python programmers as well as convenient means to enable high levels of customization.

Features
^^^^^^^^^^^^^^

   1. The data are locally stored and maintained  in a :doc:`database <database>`
   2. The database is fully accessible and modifiable through `SQLAlchemy <http://www.sqlalchemy.org>`_.  
   3. Visualization tools like heatmaps and interaction networks included
   4. The user has complete control over the information used to generate gene sets
   5. Easy to follow examples that require only a basic knowledge of Python
 
General contents:
-----------------------

.. toctree::
   :maxdepth: 1

   install
   database

Tutorials:
---------------------

.. toctree::
   :maxdepth: 1

   Database cookbook <database-cookbook>
   Using BLAST <blast>
   Gene Ontology annotation fetching <annotation-fetching>
   Gene level differential expression <deseq-example>
   Gene set analysis example <gsa-example>
   Pathway example <pathway-example>
   References <references>


Citation
^^^^^^^^^^^^

If you find ``htsint`` useful for your research or if you want to learn more about the software please refer to the following publication.

  A. J. Richards, A Herrel, & C. Bonneaud.
  `htsint: a Python library for sequencing pipelines that combines data through gene set generation <http://www.biomedcentral.com/1471-2105/16/307>`_.
  **BMC Bioinformatics**, 2015, 16, 307
