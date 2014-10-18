.. main file for htsint documentation

hts-integrate
=====================================================

``hts-integrate`` is a Python package originally conceived for the study of high-throughput sequencing data through the use of data integration approaches. However, the framework is very flexible and there are numerous other bioinformatics applications.  There is a central :doc:`database <database>` that enables :doc:`name mapping <name-mapping>`, :doc:`BLAST <blast>` querying and :doc:`annotation fetching <annotation-fetching>` via locally stored :doc:`Gene Ontology <gene-ontology>` (GO) data.  ``hts-integrate`` was created specifically to generate create custom gene-sets, based on annotation information from one or more species.  Given that these task are common to many bioinformatics pipelines, hts-integrate is not limited to sequencing pipelines and may be used in a variety of other areas of bioinformatics including, phylogenetics, gene prediction, annotation prediction, and other high-throughput analysis pipelines.

The main advantages of ``hts-integrate`` over other similar bioinformatics tools are that (1) the data is locally stored and maintained by the developer.  There is no need to rely on web-servers which can be limiting, when it comes to batch jobs, or customized queries.  (2) The database is fully accessible and modifiable by the developer through `SQLAlchemy <http://www.sqlalchemy.org>`.  This flexibility combined with a permissive license, allows forks of this project to quickly create a new database with additional tables such as PubMed documents or protein-protein interactions depending on the needs of the researchers involved.  ``hts-integrate`` also offers several levels of database access ranging from bundled convenience functions to direct SQL statements.
 
General contents:
__________________

.. toctree::
   :maxdepth: 1

   install
   database
   gene-ontology


Tutorial:
__________________

.. toctree::
   :maxdepth: 1

   blast
   name-mapping
   annotation-fetching


