#!/usr/bin/python

import time,csv,re,sys,gc
import numpy as np
from htsint.database import get_idmapping_file,get_file_sizes

timeStart = time.time()
toAdd = []
totalRecords = 0
idmappingFile = get_idmapping_file()
idmappingFid = open(idmappingFile,'rb')
reader = csv.reader(idmappingFid,delimiter="\t")


print("...gathering data this may take a few minutes")
bad = 0
current = None
noTaxa = 0
noNcbi = 0
noTaxaYesNcbi = 0

total = 0
ac2kbMap = {}
toAdd = {}

def queue_entries(toAdd):
    ## do something
    del toAdd
    gc.collect()
    toAdd = {}

for record in reader:

    continue

    #print record
    uniprotKbAc,uniprotKbEntry,ncbiId,refseq,ncbiTaxaId = None,None,None,None,None

    if len(record) != 3:
        continue

    uniprotKbAc = record[0]

    if record[1] == 'NCBI_TaxID':
        ncbiTaxaId = record[2]
    if record[1] == 'GeneID':
        ncbiId = record[2]
    if record[1] == 'UniProtKB-ID':
        uniprotKbEntry = record[2]
    if record[1] == 'RefSeq':
        refseq = record[2]

    ## update map
    if uniprotKbEntry and not ac2kbMap.has_key(uniprotKbAc):
        ac2kbMap[uniprotKbAc] = uniprotKbEntry

    ## skip the XXXX-1 like uniprot ac
    if ac2kbMap.has_key(uniprotKbAc) == False:
        continue

    ## get current key
    uniprotKbEntry = ac2kbMap[uniprotKbAc] 

    continue

    ## make new entry if necessary
    if not toAdd.has_key(uniprotKbEntry):

        ## queue in intervals
        total += 1 
        if total % 100000 == 0:
            print total
            queue_entries(toAdd)

        toAdd[uniprotKbEntry] = {'ncbi-taxa-id':None,
                                 'gene-id':None,
                                 'uniprot-ac':set([]),
                                 'refseq':set([])}
        #total += 1
        #if total == 2:
        #    print toAdd
        #    sys.exit()

    ## populate uniprot dictionary
    toAdd[uniprotKbEntry]['uniprot-ac'].update([uniprotKbAc])

    if ncbiTaxaId:
        toAdd[uniprotKbEntry]['ncbi-taxa-id'] = ncbiTaxaId
    elif ncbiId:
        toAdd[uniprotKbEntry]['gene-id'] = ncbiId
    elif refseq:
        toAdd[uniprotKbEntry]['refseq'].update([refseq])

print("run time... %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

timeStart = time.time()

#print 'toAdd loaded', len(toAdd.keys())




#print("......")
#print("total: %s"%total)
#print("noTaxa: %s"%noTaxa)
#print("noNcbi: %s"%noNcbi)
#print("noTaxaYesNcbi: %s"%noTaxaYesNcbi)
