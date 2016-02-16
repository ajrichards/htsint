#!/usr/bin/python

import time,csv,re,sys,gc
import numpy as np
from htsint.database import get_idmapping_file,get_file_sizes

timeStart = time.time()
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

ac2kbMap = {}
toAdd = {}
totalUniprot = 0#set([])

def queue_entries(toAdd,ac2kbMap):
    ## do something
    pass
    #if ac2kbMap.has_key('P07663'):
    #    print toAdd['PER_DROME']

for record in reader:

    if len(record) != 3:
        continue

    uniprotKbAc,uniprotKbEntry,ncbiId,refseq,ncbiTaxaId = None,None,None,None,None
    uniprotKbAc = record[0]

    if record[1] == 'NCBI_TaxID':
        ncbiTaxaId = record[2]
    elif record[1] == 'GeneID':
        ncbiId = record[2]
    elif record[1] == 'UniProtKB-ID':
        uniprotKbEntry = record[2]
        if not ac2kbMap.has_key(uniprotKbAc):
            ac2kbMap[uniprotKbAc] = uniprotKbEntry
    elif record[1] == 'RefSeq':
        refseq = record[2]
    else:
        continue

    ## skip the XXXX-1 like uniprot ac
    if ac2kbMap.has_key(uniprotKbAc) == False:
        continue

    ## get current key
    uniprotKbEntry = ac2kbMap[uniprotKbAc] 

    ## make new entry if necessary
    if uniprotKbEntry and not toAdd.has_key(uniprotKbEntry):

        ## queue in intervals
        totalRecords += 1 
        if totalRecords % 100000 == 0:
            queue_entries(toAdd,ac2kbMap)
            toAdd,ac2kbMap = {},{}

        #ac2kbMap[uniprotKbAc] = uniprotKbEntry
        #toAdd[uniprotKbEntry] = {'ncbi-taxa-id':None,
        #                         'gene-id':None,
        #                         'uniprot-ac':set([]),
        #                         'refseq':set([])}
        #totalUniprot.update([uniprotKbEntry])

    ## populate uniprot dictionary
    #toAdd[uniprotKbEntry]['uniprot-ac'].update([uniprotKbAc])
    #
    #if ncbiTaxaId:
    #    toAdd[uniprotKbEntry]['ncbi-taxa-id'] = ncbiTaxaId
    #elif ncbiId:
    #    toAdd[uniprotKbEntry]['gene-id'] = ncbiId
    #elif refseq:
    #    toAdd[uniprotKbEntry]['refseq'].update([refseq])

print("run time... %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

timeStart = time.time()

#print("......")
print("total: %s"%totalRecords)
#print("noTaxa: %s"%noTaxa)
#print("noNcbi: %s"%noNcbi)
#print("noTaxaYesNcbi: %s"%noTaxaYesNcbi)
