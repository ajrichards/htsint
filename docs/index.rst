.. main file for htsint documentation

hts-integrate
=====================================================

.. INCLUDE ../example_web.png 

.. figure:: example_web.png

``hts-integrate`` is a Python package used to create gene sets for the study of high-throughput sequencing data.

Features
------------

   1. The data are locally stored and maintained  in a :doc:`database <database>`
   2. The database is fully accessible and modifiable through `SQLAlchemy <http://www.sqlalchemy.org>`_.  
   3. Visualization tools like heatmaps and interaction networks included
   4. The user complete control over the information used to generate gene sets
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

   Namespace mapping <name-mapping>
   Using BLAST <blast>
   Gene Ontology annotation fetching <annotation-fetching>
   Gene level differential expression <deseq-example>
   Geneset pipeline example <pipeline-example>
