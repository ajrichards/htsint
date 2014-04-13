#!/usr/bin/env/python
"""
creates the database
"""

import sys
import numpy as np
from htsint.database import get_idmapping_file
from htsint.database import HtsDb
from htsint import __basedir__
sys.path.append(__basedir__)
try:
    from configure import CONFIG
except:
    CONFIG = None

## create the database
htsdb = HtsDb()
for dbname in htsdb.dbnames:
    if dbname in htsdb.server:
        htsdb.server.delete(dbname)
        print 'deleting', dbname
    htsdb.server.create(dbname)
htsdb.connect()

print dir(htsdb.db['htsint-uniprot'])
sys.exit()

## create a doc for each entry in the idmapping file
idmappingFile = get_idmapping_file()
idmappingFid = open(idmappingFile,'rU')
totalLines = 0
lineCount = 0
for line in idmappingFid:
    totalLines += 1
idmappingFid.close()
idmappingFid = open(idmappingFile,'rU')
allTaxa = set([])
wayPoints = [round(int(w)) for w in np.linspace(0,totalLines,100)]
toAdd = []

print 'total lines', totalLines
print "%s documents saved to 'htsint-uniprot'"%len(htsdb.db['htsint-uniprot'])

#sys.exit()

print "loading 'htsint-uniprot' database..."
for record in idmappingFid:
    record = record[:-1].split("\t")
    lineCount += 1
        
    uniprotKbAc = record[0]
    uniprotKbEntry = record[1]
    geneId = record[2]
    refseq = record[3]
    uniprotTaxon = record[13]
    doc = {'_id':uniprotKbAc,
           'uniprot-entry':uniprotKbEntry,
           'gene-id':geneId,
           'uniprot-taxa':uniprotTaxon,
           'go-terms':[]
           }

    toAdd.append(doc)
    #htsdb.db['htsint-uniprot'].save(doc,batch='ok')
    

    if lineCount in wayPoints:
        print("\t%s percent finished"%(round(lineCount/float(totalLines)*100.0)))

    if lineCount == 5:
        htsdb.db['htsint-uniprot'].bulk_save(toAdd)
        #htsdb.db['htsint-uniprot'].save(toAdd)
        break
    

idmappingFid.close()
print "%s documents saved to 'htsint-uniprot'"%len(htsdb.db['htsint-uniprot'])

#for doc in htsdb.db['htsint-uniprot']:
#    print doc, 
#    for key,item in htsdb.db['htsint-uniprot'][doc].iteritems():
#        print '\t',key,item



#print dir(db.server)
#print db.server.config()
sys.exit()


## create a connection
htsdb = HtsDb(CONFIG['dbhost'],CONFIG['dbport'])
#uniprotDb = 'htsint'

## remove the database and start fresh
dbs = htsdb.list_db()
if uniprotDb in dbs:
    htsdb.delete_db(uniprotDb)

htsdb.create_db(uniprotDb)

## read in the files
#geneInfo = read_gene_info_file()

## put all of the genes from the mapping files into the db
idmappingFile = get_idmapping_file()
idmappingFid = open(idmappingFile,'rU')
result = {}
allTaxa = set([])

debug = 0
for record in idmappingFid:
    record = record[:-1].split("\t")
    debug += 1
        
    uniprotKbAc = record[0]
    uniprotKbEntry = record[1]
    geneId = record[2]
    refseq = record[3]
    uniprotTaxon = record[13]
    #result[uniprotKbAc] = [geneId,uniprotKbEntry,refseq]

    doc = """
    {
        "value":
        {
            "UniProtKbEntry":"%s",
            "GeneID":"%s",
            "RefSeq":"%s",
            "Tags":["plankton", "baseball", "decisions"],
            "Body":"I decided today that I don't like baseball. I like plankton."
        }
    }
    """%(uniprotKbEntry,geneId,refseq)
    htsdb.save_doc(uniprotDb, doc, uniprotKbAc)

    if debug == 2:
        break

idmappingFid.close()

g = htsdb.open_doc(uniprotDb,uniprotKbAc)
print dir(g)
print type(g)
print g
#htsdb.list_doc(uniprotDb)
#htsdb.info_db(uniprotDb)


#delete_db(dbName)




'''
    print "\nCreate a document 'mydoc' in database 'mydb':"
    doc = """
    {
        "value":
        {
            "Subject":"I like Planktion",
            "Author":"Rusty",
            "PostedDate":"2006-08-15T17:30:12-04:00",
            "Tags":["plankton", "baseball", "decisions"],
            "Body":"I decided today that I don't like baseball. I like plankton."
        }
    }
    """
    foo.saveDoc('mydb', doc, 'mydoc')

    print "\nCreate a document, using an assigned docId:"
    foo.saveDoc('mydb', doc)

    print "\nList all documents in database 'mydb'"
    foo.listDoc('mydb')

    print "\nRetrieve document 'mydoc' in database 'mydb':"
    foo.openDoc('mydb', 'mydoc')

    #print "\nDelete document 'mydoc' in database 'mydb':"
    #foo.deleteDoc('mydb', 'mydoc')

    print "\nList all documents in database 'mydb'"
    foo.listDoc('mydb')

    print "\nList info about database 'mydb':"
    foo.infoDb('mydb')

    print "\nDelete database 'mydb':"
    foo.deleteDb('mydb')

    print "\nList databases on server:"
    foo.listDb()

if __name__ == "__main__":
    test()
'''
