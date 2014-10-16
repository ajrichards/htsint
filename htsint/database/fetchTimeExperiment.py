#!/usr/bin/python

import sys,time
from sqlalchemy.sql import select
from htsint.database import db_connect
from htsint.database import Taxon,taxa_mapper,Gene,gene_mapper 
session,engine = db_connect()
conn = engine.connect()

widget = Gene#Taxon
print("scanning %s"%widget.__tablename__)

timeStart = time.time()
myDict = {}
s = select([widget.id,widget.ncbi_id])
_result = conn.execute(s)  
result = [row for row in _result]
print("core: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
print result[:5]
print result[0]['ncbi_id']
sys.exit()

timeStart = time.time()
for t in session.query(widget).yield_per(5):
    myDict[t.ncbi_id] = t.id
print("yield per: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
