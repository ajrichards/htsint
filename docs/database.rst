.. main file for lpedit documentation

Database
======================

Because ``htsint`` was originally intended for use with non-model organisms and because it is useful to make inference across multiple species the database is constructed around the information contained within the `UniProt <http://www.uniprot.org>`_ Knowledgebase (UniProtKB).  For each sequence included in UniProtKB an *accession number* is assigned.  These numbers are consistent across releases.  There is also an *Entry Name* or *ID* that is used by UniProtKB.  An *Entry Name* is usually biologically meaningful, however, it may change over time and multiple accession numbers can be associated with a single *Entry Name*.  UniProt has a `page with more infomation <http://www.uniprot.org/faq/6>`_.   

Because *Entry Name* is not a stable identifier, we begin our database with UniProtKB *accession numbers*.  We are also interested in :doc:`name mapping <name-mapping>` so we include a table for `GeneIDs <http://www.ncbi.nlm.nih.gov/gene>`_.  UniProtKB accessions and NCBI's GeneID have an associated taxon.  Because the integration of functional annotation information is crucial to genomics applications we include tables for `Gene Ontology <http://www.geneontology.org>`_ terms and annotations.


Database migration or to create a backup
----------------------------------------

   .. code-block:: bash

      ~$ pg_dump -h localhost -U dbuser dbname > htsint.sql

