#!/usr/bin/python




from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene,Uniprot,Refseq

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


results = conn.execute(Uniprot.__table__.select(Uniprot.uniprot_entry.in_(uniprotEntries)))
upEntry2Gene = dict([(str(row.uniprot_entry),str(row.gene_id)) for row  in results])
upEntry2Taxa = dict([(str(row.uniprot_entry),str(row.taxa_id)) for row  in results])

for key, item in upEntry2Gene.iteritems():
    print key, item


#def row2dict(row):
#    d = {}
#    for column in row.__table__.columns:
#        d[column.name] = str(getattr(row, column.name))

#    return d


#query = session.query(Uniprot).filter_by(uniprot_ac='A8WGU8').all()
#print query
