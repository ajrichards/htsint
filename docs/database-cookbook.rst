.. main file for lpedit documentation

Common database tasks
========================

This pages serves a a list of useful examples demonstrating how to interact with the :doc:`database`.  First, we need to establish a connection with the database.

>>> from htsint.database import db_connect()
>>> session, engine = db_connect()

Fetching common names associated with scientific names
---------------------------------------------------------

>>> from htsint.database import Taxon
>>> taxaList = ['Thamnophis sirtalis',\
                'Hemiaspis signata',\
                'Pantherophis guttatus',\
                'Crotalus horridus']
>>> tQuery = session.query(Taxon).filter(Taxon.name.in_(taxaList)).all()

>>> for tq in tQuery:
    print tq.ncbi_id, tq.name, tq.common_name_1
8665 Ophiophagus hannah king cobra
8729 Crotalus adamanteus eastern diamondback rattlesnake
35024 Crotalus horridus timber rattlesnake
8637 Micrurus fulvius eastern coral snake

Get number of protein coding genes
--------------------------------------

>>> from htsint.database import Uniprot
>>> tqueryId = session.query(Taxon).filter(Taxon.ncbi_id=='7227').first().id
>>> uniprotQuery = session.query(Uniprot).filter_by(taxa_id=tqueryId).all()
>>> _codingGenes = list(set([u.gene_id for u in uniprotQuery]))
>>> codingGenes = [g.ncbi_id for g in session.query(Gene).filter(Gene.id.in_(_codingGenes)).all()]
>>> print("coding genes: %s"%(len(codingGenes)))
coding genes: 13587

Load gene information as a dictionary
------------------------------------------

For large queries, the fastest way to access the database is usually using the SQlAlchemy core syntax.

>>> from htsint.database import Gene
>>> from sqlalchemy.sql import select
>>> conn = engine.connect()

>>> tqueryId = session.query(Taxon).filter(Taxon.ncbi_id=='7227').first().id
>>> s = select([Gene.ncbi_id,Gene.description,Gene.symbol],Gene.taxa_id==tqueryId)
>>> _geneQueries = conn.execute(s)
>>> geneQueries = _geneQueries.fetchall()
>>> gene2desc = dict([(str(r['ncbi_id']),str(r['description'])) for r in geneQueries])
>>> gene2sym = dict([(str(r['ncbi_id']),str(r['symbol'])) for r in geneQueries])

>>> geneIds = ['44817','10216631','34790']
>>> for gid in geneIds:
        print("gene: %s (%s) - %s"%(gene2sym[gid],gid,gene2desc[gid]))
gene: for (44817) - foraging
gene: caboose (10216631) - caboose
gene: Sos (34790) - Son of sevenless
