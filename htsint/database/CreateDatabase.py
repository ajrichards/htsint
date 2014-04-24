#!/usr/bin/env python
"""
creates htsint database
(1) Edit/create the configure.py file and specify the database parameters
(2) Run 'FetchDbData.py'
(3) Run this file.
A logfile of the output is automatically created in the 'data' directory.
Use TestDatabase.py to make sure tables were filled correctly.
There is also a unittest for the database and the examples within show some
of the different ways of interacting with the database.
For more information about interacting with databases using SQLalchemy
    http://docs.sqlalchemy.org/en/rel_0_8/orm/tutorial.html
Because gene names are added/removed/changed on a regular basis
it is best to re-run FetchDbdata.py and then this script.
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
from DatabaseTools import populate_go_terms, populate_go_annotations,get_taxa_list
from GeneOntologyLib import read_annotation_file,get_annotation_file

## prepare a log file
fid = open(os.path.join(CONFIG['data'],'createdb.log'),'wa')
writer = csv.writer(fid)

def push_out(line):
    writer.writerow([line])
    print(line)

push_out(sys.argv[0])
push_out(time.asctime())
push_out("Getting ready to create database...")

## conect to the database
session,engine = db_connect(verbose=False)
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)

push_out("Creating database with...")
for t in Base.metadata.sorted_tables:
   push_out("\t"+t.name)

## get a list of geneids from uniprot
timeStart = time.time()
push_out("extracting gene ids...")
geneIds, idmapLineCount = get_geneids_from_idmapping()
push_out('%s geneIds were found in the idmapping file'%len(geneIds))
sys.exit()
push_out("extracting taxa list...")
taxaList = get_taxa_list()
print(len(taxaList))
push_out("...total time taken to extract data: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart)))



## taxa table
push_out("Populating the database with %s taxa"%len(taxaList))
timeStr,addedStr = populate_taxon_table(taxaList,engine)
push_out(timeStr)
push_out(addedStr)

## gene table
push_out("Populating the database with %s genes"%len(geneIds.keys()))
timeStr,addedStr = populate_gene_table(geneIds,taxaList,session,engine)
push_out(timeStr)
push_out(addedStr)

##  uniprot table
push_out("Populating the database with %s uniprot entries"%(idmapLineCount))
timeStr,addedStr = populate_uniprot_table(idmapLineCount,geneIds,session,engine)
push_out(timeStr)
push_out(addedStr)

## populate the go-terms
push_out("Populating the database with for go terms...")
timeStr,addedStr = populate_go_terms(engine)
push_out(timeStr)
push_out(addedStr)

## populate the go-annotations
push_out("Populating the database with for go annotations...")
timeStr,addedStr = populate_go_annotations(session,engine)
push_out(timeStr)
push_out(addedStr)


print_db_summary()
fid.close()
