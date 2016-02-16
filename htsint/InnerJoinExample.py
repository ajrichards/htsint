import time,sys
from htsint.database import db_connect,Gene,GoAnnotation,GoTerm
from htsint.database import fetch_annotations,gene_mapper


## variables
session,engine = db_connect()
geneList = ['30970']#,'30971','30972','30973','30975']
expEvidCodes = ["EXP","IDA","IPI","IMP","IGI","IEP"]
compEvidCodes = ["ISS","ISO","ISA","ISM","IGC","RCA"]
statEvidCodes = ["TAS","NAS","IC"]
nonCuratedEvidCodes = ["IEA"]
acceptedCodes = expEvidCodes + statEvidCodes
annotations = {}
aspect = 'biological_process'

timeStart = time.time()
geneQueries = session.query(Gene).filter(Gene.ncbi_id.in_(geneList)).all()
print("...extraction time 1: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

#timeStart = time.time()
#geneIdMap = gene_mapper(session,ncbiIdList=geneList)
#print("...extraction time 1: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

## get all results
timeStart = time.time()
for geneQuery in geneQueries:
    annotations[geneQuery.ncbi_id] = set([])
    print geneQuery.ncbi_id
    annotations[geneQuery.ncbi_id].update(session.query(GoAnnotation).filter_by(gene_id=geneQuery.id).all())
print("...extraction q1: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

for key,items in annotations.iteritems():
    annotations[key] = list(items)
    if None in items:
        annotations[key].remove(None)

for at in annotations['30970']:
    print at,session.query(GoTerm).filter_by(id = at.go_term_id).first().aspect

## get results with table join by aspect
timeStart = time.time()

print 'filtered results'
for geneQuery in geneQueries:
    results = session.query(GoAnnotation).join(GoTerm).\
              filter(GoAnnotation.gene_id==geneQuery.id).\
              filter(GoAnnotation.evidence_code.in_(['ISS'])).\
              filter(GoTerm.aspect=='molecular_function').all()
    
    for r in results:
        print r, session.query(GoTerm).filter_by(id = r.go_term_id).first().aspect
