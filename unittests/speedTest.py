#!/usr/bin/python



import time
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene,Uniprot,Refseq
from htsint.database import uniprot_mapper

session,engine = db_connect()
conn = engine.connect()

uniprotEntries = ["KCNQ4_MOUSE","CSMT1_XENTR","CSMT1_MOUSE","MILK2_MOUSE","MILK2_RAT",
                  "MILK1_RAT","MILK1_MOUSE","MICA3_MOUSE","MCA3A_DANRE","MCA3B_DANRE",
                  "MICA1_DANRE","MCA2B_DANRE","MICLK_MOUSE","MICLK_RAT","MILK2_RAT",
                  "MILK2_MOUSE","MILK1_RAT","MILK1_MOUSE","MICA3_MOUSE","MICA2_RAT",
                  "MCA3A_DANRE","MICA2_MOUSE","EHBP1_MOUSE","EH1L1_MOUSE","SPTB2_MOUSE",
                  "MCA2B_DANRE","SPTN2_RAT","MICA2_XENTR","SPTCB_DROME","ACTN_DROME",
                  "SPTB1_MOUSE","ACTN2_MOUSE","ACTN3_MOUSE","ACTN2_CHICK","ACTN1_RAT",
                  "ACTN1_CHICK","ACTN1_MOUSE","CYTSA_CHICK","MCA3B_DANRE","CYTSA_CANFA",
                  "CYTSA_DANRE","CYTSA_MOUSE","AIN1_SCHPO","MICA1_DANRE","CYTSA_RAT",
                  "SYNE2_MOUSE","ACTN4_CHICK","ACTN4_MOUSE","ACTN4_RAT","CYTSA_XENTR",
                  "CYTSB_MOUSE","SMTL2_MOUSE","SMTN_MOUSE","DYST_MOUSE","PLEC_RAT",
                  "PLEC_MOUSE","DMD_CHICK","DMD_CANFA","DMD_MOUSE","MICA1_RAT",
                  "SMTL1_MOUSE","MICA1_MOUSE","MACF1_MOUSE","MACF1_RAT","DMD_CAEEL",
                  "MILK2_MOUSE","MILK2_RAT","MILK1_RAT","MILK1_MOUSE","ACTN4_CHICK"
                  "ACTN4_RAT","ACTN4_MOUSE","ACTN1_CHICK","ACTN_DROME","ACTN1_RAT"
                  "ACTN1_MOUSE","ACTN3_MOUSE","SPTCB_DROME","ACTN2_MOUSE","ACTN2_CHICK"]

## using select method
timeStart = time.time()
s = select([Uniprot.uniprot_entry,Uniprot.taxa_id,Uniprot.gene_id]).where(Uniprot.uniprot_entry.in_(uniprotEntries))
_upQueries = conn.execute(s)
upQueries = _upQueries.fetchall()
upEntry2Gene = dict([(str(uquery['uniprot_entry']),str(uquery['gene_id'])) for uquery in upQueries])
upEntry2Taxa = dict([(str(uquery['uniprot_entry']),str(uquery['taxa_id'])) for uquery in upQueries])
print("Test 1: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## using direct table access method
timeStart = time.time()
results = conn.execute(Uniprot.__table__.select(Uniprot.uniprot_entry.in_(uniprotEntries)))
upEntry2Gene, upEntry2Taxa = {},{}
for row in results:
    upEntry2Gene[str(row.uniprot_entry)] = str(row.gene_id)
    upEntry2Taxa[str(row.uniprot_entry)] = str(row.taxa_id)

for key, item in upEntry2Gene.iteritems():
    print key, item
print("Test 2: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## using htsint's mapper
timeStart = time.time()
uMapper = uniprot_mapper(session,uniprotIdList=uniprotEntries,gene=True,taxa=True)
print("Test 3: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
print(len(uMapper.keys()))
print(uMapper[uniprotEntries[0]])
