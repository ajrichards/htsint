#!/usr/bin/python

import sys,time
from htsint.database import db_connect
from htsint.database import Taxon,taxa_mapper,Gene,gene_mapper 
session,engine = db_connect()

from sqlalchemy import and_, func

def column_windows(session, column, windowsize):
    """Return a series of WHERE clauses against 
    a given column that break it into windows.

    Result is an iterable of tuples, consisting of
    ((start, end), whereclause), where (start, end) are the ids.

    Requires a database that supports window functions, 
    i.e. Postgresql, SQL Server, Oracle.

    Enhance this yourself !  Add a "where" argument
    so that windows of just a subset of rows can
    be computed.

    """
    def int_for_range(start_id, end_id):
        if end_id:
            return and_(
                column>=start_id,
                column<end_id
            )
        else:
            return column>=start_id

    q = session.query(
                column, 
                func.row_number().\
                        over(order_by=column).\
                        label('rownum')
                ).\
                from_self(column)
    if windowsize > 1:
        q = q.filter("rownum %% %d=1" % windowsize)

    intervals = [id for id, in q]

    while intervals:
        start = intervals.pop(0)
        if intervals:
            end = intervals[0]
        else:
            end = None
        yield int_for_range(start, end)

def windowed_query(q, column, windowsize):
    """"Break a Query into windows on a given column."""

    for whereclause in column_windows(q.session, 
                                      column, windowsize):
        for row in q.filter(whereclause).order_by(column):
            yield row



widget = Gene

print("scanning %s"%widget.__tablename__)

## yeild per
timeStart = time.time()
for t in session.query(widget).yield_per(5):
    pass
    #if ncbiIdList and not ncbiIdList.has_key(str(t)):
    #    continue
print("Yield per: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

timeStart = time.time()
q = session.query(widget)
myDict = {}
for t in windowed_query(q,widget.ncbi_id,1000):
    pass
    #myDict[str(t.ncbi_id)] = t.id
print("Windowed query: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))





sys.exit()
#q = s.query(Widget)
#
#for widget in windowed_query(q, Widget.data, 1000):
#    print "data:", widget.data



def page_query(q):
    offset = 0
    while True:
        r = False
        for elem in q.limit(1000).offset(offset):
           r = True
           yield elem
        offset += 1000
        if not r:
            break

"""
## taxa
timeStart = time.time()
myDict = {}
for t in page_query(session.query(Taxon)):
    myDict[str(t.ncbi_id)] = t.id
print("Page query: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## yield per 5
timeStart = time.time()
taxaIdMap = taxa_mapper(session,yieldPer=5)
print("Yield per 5: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
"""

## genes
timeStart = time.time()
myDict = {}
for t in page_query(session.query(Gene)):
    myDict[str(t.ncbi_id)] = t.id
print("Page query: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## yield per 5
#timeStart = time.time()
#taxaIdMap = gene_mapper(session,yieldPer=5)
#print("Yield per 5: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))







sys.exit()
#ncbiIdList = taxaIdMap.keys()[1000]
#jack = session.query(Taxa).\
#       options(joinedload('addresses')).\
#       filter_by(name='jack').all() 


#print taxaIdMap.


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
