#!/usr/bin/python

import time,csv,re,sys
import numpy as np
from htsint.database import get_idmapping_file,get_file_sizes




idmapCount,geneInfoCount = get_file_sizes()
print 'idmapcount', idmapCount
print 'gene info count', geneInfoCount

sys.exit()


timeStart = time.time()
toAdd = []
totalRecords = 0
idmappingFile = "/usr/local/share/htsint/idmapping.dat.db"
idmappingFid = open(idmappingFile,'rb')
reader = csv.reader(idmappingFid,delimiter="\t")

print("...gathering data this may take a few minutes")
bad = 0
current = None
uniprotKbEntry,ncbiId,refseq,ncbiTaxaId = None,None,None,None

noTaxa = 0
noNcbi = 0
total = 0
noTaxaYesNcbi = 0
weirdResults = 0

for record in reader:

    if len(record) != 3:
        continue

    uniprotKbAc = record[0]

    if current == None:
        current = uniprotKbAc

    if record[1] == 'NCBI_TaxID':
        ncbiTaxaId = record[2]
    if record[1] == 'GeneID':
        ncbiId = record[2]
    if record[1] == 'UniProtKB-ID':
        uniprotKbEntry = record[2]
    if record[1] == 'RefSeq':
        refseq = record[2]

    if re.search("\t",record[2]):
        print record
        weirdResults += 1

    ## check to see if entry is finished
    if current != uniprotKbAc:
        current = uniprotKbAc
        total += 1

        if ncbiTaxaId == None:
            noTaxa += 1
        if ncbiId == None:
            noNcbi += 1
        if ncbiId and ncbiTaxaId == None:
            noTaxaYesNcbi += 1

    
        #if total > 100:
        #    sys.exit()

            

        #print "\nuniprotid:%s\nncbiid%s\nrefseq:%s\ntaxaid:%s"%(uniprotKbEntry,ncbiId,refseq,ncbiTaxaId)

        #if None in [uniprotKbEntry,ncbiId,refseq,ncbiTaxaId]:
        #    print "\nuniprotid:%s\nncbiid%s\nrefseq:%s\ntaxaid:%s"%(uniprotKbEntry,ncbiId,refseq,ncbiTaxaId)
        #else:
        #    print "...", uniprotKbEntry
        #sys.exit()

        uniprotKbEntry,ncbiId,refseq,ncbiTaxaId = None,None,None,None


print("......")
print("total: %s"%total)
print("noTaxa: %s"%noTaxa)
print("noNcbi: %s"%noNcbi)
print("noTaxaYesNcbi: %s"%noTaxaYesNcbi)
