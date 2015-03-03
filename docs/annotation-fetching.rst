.. main file for lpedit documentation

Annotation fetching
=======================

A very common bioinformatics task is that of obtaining current annotations for a given gene or protein.  There is a great deal of flexibility in how annotations may be obtained from the database and some strategies are more or less appropriate for a specific application.  First, however it is important that you understand how the `Gene Ontology (GO) <http://geneontology.org>`_ information is stored in the :doc:`database <database>`.


Fetching annotations with a list of NCBI or UniProt identifiers
-----------------------------------------------------------------

First, we need to connect to the database.

>>> from htsint.database import fetch_annotations,db_connect
>>> session, engine = db_connect()
>>> conn = engine.connect()

The ``conn`` or connection is for usage with the `SQLAlchemy core API <http://docs.sqlalchemy.org/en/rel_0_9/core/>`_ and is not necessary for basic annotation fetching.  If you have a list of NCBI gene identifiers or UniProt identifiers (UniprotAc) then the easiest way to get the annotations is to

>>> annotations = fetch_annotations(['31251'],engine,idType='ncbi',useIea=False,aspect='biological_process')
>>> print(annotations)
>>> print(annotations['31251'][0])
('GO:0042752', 'regulation of circadian rhythm')

Simply use `idType='uniprot'` if the ids are from that namespace.  As you can see the results are contained in a dictionary with a key for each NCBI identifier.  The results associated with each identifier are a unique list of tuples each containing the GO identifier and GO name. So if you wanted only the GO names simply use a list comprehension.

   .. code-block:: python

      print([a[1] for a in annotations['31251']])

A results list can be quickly sorted with:

   .. code-block:: python

      sortby_inplace(annotations['31251'], 1)

.. note::
  For large lists it may be faster to search for all annotations using ``fetch_taxa_annotations`` (see below) then to trim the list


Fetching directly with SQLAlchemy
-------------------------------------

There are two API's from SQLAlchemy the ORM and the Core.

Using the ORM
^^^^^^^^^^^^^^^^^^^^

   >>> subq1 = session.query(GoTerm).\
               filter(GoTerm.aspect==aspect).\
               subquery()
   >>> subq2 = session.query(GoAnnotation).\
             filter(GoAnnotation.gene_id==gquery['id']).\
             subquery()

   >>> q = session.query(subq1.c.name).\
           filter(subq1.c.id==subq2.c.go_term_id).all()
   

Using the Core
^^^^^^^^^^^^^^^^^^
     
   >>> s = select([GoTerm.name],GoAnnotation.go_term_id==GoTerm.id).\
         where(GoAnnotation.gene_id == gquery['id']).\
         where(GoTerm.aspect==aspect)

   >>> _results = conn.execute(s)
   >>>  results = [tuple(map(str,r)) for r in _results.fetchall()]   

Fetching annotations by taxa
---------------------------------------------------------

The syntax is similar to individual gene/protein identifiers.

>>> from htsint.database import fetch_taxa_annotations,db_connect
>>> session, engine = db_connect()
>>> conn = engine.connect()
>>> geneAnnots,uniprotAnnots = fetch_taxa_annotations(['7091'],engine,useIea=False,verbose=True)
>>> print(uniprotAnnots['Q9NL89'])
['GO:0045088']
