#!/usr/bin/env python
"""
Creates the basic shell of the database
Database is populated afterwards with PopulateDatabase.py

We populate the database with a several model organisms
including:

 * Drosophila melanogaster  (7227)
 * Mus musculus             (10090)
 * Homo sapiens             (9606)
 * Saccharomyces cerevisiae (4932)

A logfile of the output is automatically created

Next:
    Use TestDatabase.py to make sure tables were filled correctly

For more information about interacting with databases using SQLalchemy
    http://docs.sqlalchemy.org/en/rel_0_8/orm/tutorial.html

Because gene names are added/removed/changed on a regular basis
It is easier when updating to re-run this script and then follow up
with the PopulateDatabase.py afterwards. 
"""

### make imports
import sys,os,re,time,csv
from htsint import __basedir__

sys.path.append(__basedir__)
try:
    from configure import CONFIG
except:
    CONFIG = None

from DatabaseTables import Base,Taxon,Gene,Uniprot,GoTerm,GoAnnotation
from DatabaseTools import db_connect, get_geneids_from_idmapping,print_db_summary
from DatabaseTools import populate_taxon_table,populate_gene_table,populate_uniprot_table
from DatabaseTools import populate_go_terms, populate_go_annotations
from GeneOntologyLib import read_annotation_file,get_annotation_file

## prepare a log file
fid = open(os.path.join(CONFIG['data'],'createdb.log'),'w')
writer = csv.writer(fid)

def push_out(line):
    writer.writerow([line])
    print(line)

push_out(sys.argv[0])
push_out(time.asctime())
push_out("Getting ready to create database...")

## conect to the database
session,engine = db_connect(verbose=False)
#print dir(Base.metadata)
#print Base.metadata.tables.keys()
#print dir(Base.metadata.tables)
#sys.exit()
#print Base.metadata.sorted_tables
#Uniprot.__table__.drop(engine)

session.commit()

#Base.metadata.remove(GoAnnotation)
#print dir(Base.metadata)
#Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)

push_out("Creating database with...")
for t in Base.metadata.sorted_tables:
   push_out("\t"+t.name)

## get a list of geneids from uniprot
geneIds, idmapLineCount = get_geneids_from_idmapping()
push_out('%s geneIds were found in the idmapping file'%len(geneIds))

## use the list of geneids and the annotation file to to get the taxa list
geneInfoFile = os.path.join(CONFIG['data'],"gene_info.db")
geneInfoFid = open(geneInfoFile,'rU')
taxaList = set([])
header = geneInfoFid.next()
for record in geneInfoFid:
    record = record.rstrip("\n")
    record = record.split("\t")
    if re.search("^\#",record[0]):
        continue
    taxaList.update([record[0]])

push_out('%s elements extracted from gene info file'%len(geneIds.keys()))
annotationFile = get_annotation_file()
annotationFid = open(annotationFile,'rU')
annotsCount = 0
annotatedIds = {}

for record in annotationFid:
    record = record[:-1].split("\t")
    if record[0][0] == "!":
        continue
    if record[0] != 'UniProtKB':
        continue

    annotsCount += 1
    taxon = re.sub("taxon:","",record[12])
    if taxon == "" or re.search("\|",taxon):
        continue

    annotatedIds[record[1]] = None
    taxaList.update([taxon])

taxaList = list(taxaList)

## taxa table
#push_out("Populating the database with %s taxa"%len(taxaList))
#timeStr,addedStr = populate_taxon_table(taxaList,session,engine)
#push_out(timeStr)
#push_out(addedStr)

## gene table
#push_out("Populating the database with %s genes"%len(geneIds.keys()))
#timeStr,addedStr = populate_gene_table(geneIds,session)
#push_out(timeStr)
#push_out(addedStr)

##  uniprot table
push_out("Populating the database with %s uniprot entries"%(idmapLineCount))
timeStr,addedStr = populate_uniprot_table(idmapLineCount,session)
push_out(timeStr)
push_out(addedStr)

## populate the go-terms
#push_out("Populating the database with for go terms...)
#timeStr,addedStr = populate_go_terms(session)
#push_out(timeStr)
#push_out(addedStr)

## populate the go-annotations
#push_out("Populating the database with for go terms...)
#timeStr,addedStr = populate_go_terms(session)
#push_out(timeStr)
#push_out(addedStr)



## go terms and go annotations tables
#push_out("Populating the database with %s annotations"%(annotCount))
#timeStr,addedStr = populate_go_tables(session)
#push_out(timeStr)
#push_out(addedStr)


print_db_summary()
fid.close()
