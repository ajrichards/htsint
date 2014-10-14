#!/usr/bin/python

import sys,time
from htsint.database import db_connect
from htsint.database import taxa_mapper 
session,engine = db_connect()


## yield per 5
timeStart = time.time()
taxaIdMap = taxa_mapper(session,yieldPer=5)
print("Yield per 5: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## yield per 500
#timeStart = time.time()
#taxaIdMap = taxa_mapper(session,yieldPer=500)
#print("Yield per 500: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## yield per 5000
#timeStart = time.time()
#taxaIdMap = taxa_mapper(session,yieldPer=5000)
#print("Yield per 5000: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## yield per 5 with subset
timeStart = time.time()
ncbiIdList = taxaIdMap.keys()[5000]
taxaIdMap = taxa_mapper(session,ncbiIdList=ncbiIdList,yieldPer=5)
print("Yield per 5 with subset: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))




sys.exit()

#goid GO:0005737
#evidencecode IDA
#pubmed PMID:11168596
#uniprotid Q9GNQ1
#geneid None
#taxon 146127


tQuery = session.query(Taxon).filter_by(ncbi_id='146127').first()
uQuery = session.query(Uniprot).filter_by(uniprot_ac='Q9GNQ1').first()
gQuery = session.query(Gene).filter_by(id=uQuery.gene_id).first()

print 'tq',tQuery
print 'uq',uQuery

print 'ac',uQuery.uniprot_ac
print 'entry',uQuery.uniprot_entry
print 'refseq', uQuery.refseq
print 'taxa_id', uQuery.taxa_id
print 'gene_id', uQuery.gene_id

print 'gene', gQuery.ncbi_id


tax = session.query(Taxon).filter_by(id=uQuery.taxa_id).first()
print tax.ncbi_id, tax.name, tax.common_name_1

sys.exit()













print dir(taxaIdMap)
print len(taxaIdMap.keys())
if 2357534 in taxaIdMap.values():
    print True
else:
    print False
