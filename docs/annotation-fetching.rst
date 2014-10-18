.. main file for lpedit documentation

Annotation fetching
=======================

A very common bioinformatics task is that of obtaining current annotations for a given gene or protein.  There is a great deal of flexibility in how annotations may be obtained from the database and some strategies are more or less appropriate for a specific application.  First, however it is important that you understand how the :doc:`Gene Ontology <gene-ontology>` information is stored in the :doc:`database <database>`.


Fetching annotations with a list of NCBI or UniProt identifiers
______________________________________________________________________


First, we need to connect to the database.

   > from htsint.database import fetch_annotations,db_connect
   > session, engine = db_connect()
   > conn = engine.connect()

The ``conn`` or connection is for usage with the `SQLAlchemy core API <http://docs.sqlalchemy.org/en/rel_0_9/core/>`_ and is not necessary for basic annotation fetching.  If you have a list of NCBI gene identifiers or UniProt identifiers (UniprotAc) then the easiest way to get the annotations is to

   .. code-block:: python

      from htsint.database import fetch_annotations

      annotations = fetch_annotations(['31251'],engine,idType='ncbi',useIea=False,aspect='biological_process')
      print(annotations)

Simply use `idType='uniprot'` if the ids are from that namespace.  As you can see the results are contained in a dictionary with a key for each NCBI identifier.  The results associated with each identifier are a unique list of tuples each containing the GO identifier and GO name. So if you wanted only the GO names simply use a list comprehension.

   > print([a[1] for a in annotations['31251']])

A results list can be quickly sorted say by name using:

   > sortby_inplace(annotations['P56645'], 1)


Fetching directly with SQLAlchemy
______________________________________


There are two API's from SQLAlchemy the ORM and the Core.

Using the ORM
^^^^^^^^^^^^^^^^^^^^

   > subq1 = session.query(GoTerm).\
             filter(GoTerm.aspect==aspect).\
             subquery()
   > subq2 = session.query(GoAnnotation).\
             filter(GoAnnotation.gene_id==gquery['id']).\
             subquery()

   > q = session.query(subq1.c.name).\
         filter(subq1.c.id==subq2.c.go_term_id).all()
   

Using the Core
^^^^^^^^^^^^^^^^^^
     
   > s = select([GoTerm.name],GoAnnotation.go_term_id==GoTerm.id).\
         where(GoAnnotation.gene_id == gquery['id']).\
         where(GoTerm.aspect==aspect)

   > _results = conn.execute(s)
   >  results = [tuple(map(str,r)) for r in _results.fetchall()]   

