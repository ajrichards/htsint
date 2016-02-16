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


Load more gene and GO information that you will likely need
--------------------------------------------------------------

If your query is not specific and you will need to make many genome wide queries it can be helpful to load all gene and GO information into memory at once.

.. code-block:: python

    from htsint.database import fetch_taxa_annotations,db_connect
    from htsint.database import Gene,Taxon,GoTerm

    ## connect to the db
    session,engine = db_connect()
    conn = engine.connect()

    ## load gene centric info
    tqueryId = session.query(Taxon).filter(Taxon.ncbi_id=='9606').first().id
    s = select([Gene.ncbi_id,Gene.description,Gene.symbol,Gene.map_location,Gene.chromosome],Gene.taxa_id==tqueryId)
    _geneQueries = conn.execute(s)
    geneQueries = _geneQueries.fetchall()
    gene2desc = dict([(str(r['ncbi_id']),str(r['description'])) for r in geneQueries])
    sym2gene = dict([(str(r['symbol']),str(r['ncbi_id'])) for r in geneQueries])
    gene2maploc =  dict([(str(r['ncbi_id']),str(r['map_location'])) for r in geneQueries])
    gene2chrom =  dict([(str(r['ncbi_id']),str(r['chromosome'])) for r in geneQueries])

    ## load go term defs
    aspectList = ['biological_process','molecular_function']
    _termDicts = dict([(aspect, session.query(GoTerm).filter_by(aspect=aspect).all()) for aspect in aspectList])
    term2name,term2desc = {},{}

    for key,terms in _termDicts.iteritems():
        term2name[key] = dict([(term.go_id,term.name) for term in terms])
        term2desc[key] = dict([(term.go_id,term.description) for term in terms])

    print(term2name['biological_process']['GO:1900308'])
    print(term2desc['biological_process']['GO:1900308'])
		
    geneAnnotsBP,uniprotAnnotsBP = fetch_taxa_annotations(['9606'],engine,useIea=True,aspect='biological_process')
    geneAnnotsMF,uniprotAnnotsMF = fetch_taxa_annotations(['9606'],engine,useIea=True,aspect='molecular_function')
