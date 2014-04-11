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
from DatabaseTools import db_connect, read_gene_info_file
from DatabaseTools import populate_taxon_table,populate_gene_table,populate_uniprot_table,populate_go_tables
from GeneOntologyLib import read_idmapping_file,read_annotation_file

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
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)

push_out("Creating database with...")
for t in Base.metadata.sorted_tables:
   push_out("\t"+t.name)

## read the mappings and annotations
print("reading necessary files into memory")
annotations = read_annotation_file()
mappings = read_idmapping_file()
geneInfo = read_gene_info_file()
print("loaded")

push_out('Total annotations = %s'%(len(annotations.keys())))
push_out('Total mappings = %s'%(len(mappings.keys())))

## get genes and taxa in mappings
taxaList = set([])

## replace the multiple geneIds in mappings with corresponding geneIds from gene_info
for uniprotac, uniprotmap in mappings.iteritems():
    geneId = uniprotmap[0]
    
    ## use only a gene id that matches
    if re.search(";",geneId):
        for _geneId in [re.sub("\s+","",gid) for gid in geneId.split(";")]:
            if geneInfo.has_key(_geneId):
                geneId = _geneId
                mappings[uniprotac][0] = geneId

geneIdList = [g[0] for g in mappings.values()]
for geneId in geneIdList:
    if geneId == "":
        continue
    
    ## use only the first match
    if not geneInfo.has_key(geneId):
        print("WARNING cannot find %s in gene_info file"%geneId)
        continue

    taxaList.update([geneInfo[geneId][0]])

taxaList = list(taxaList)

## add all taxa present in the mappings file
push_out("Populating the database with %s taxa"%len(taxaList))
timeStr,addedStr = populate_taxon_table(taxaList,session)
push_out(timeStr)
push_out(addedStr)

## gene table
push_out("Populating the database with %s genes"%len(mappings.keys()))
timeStr,addedStr = populate_gene_table(mappings,geneInfo,session)
push_out(timeStr)
push_out(addedStr)

##  table
push_out("Populating the database with %s uniprot entries"%len(mappings.keys()))
timeStr,addedStr = populate_uniprot_table(mappings,annotations,session)
push_out(timeStr)
push_out(addedStr)


sys.exit()

## taxon table
goTaxa = get_all_go_taxa()
taxaList = ["7227"] + CONFIG['taxa']
taxaList = taxaList + goTaxa
push_out("Attempting to populate the database with %s taxa"%(len(taxaList)))

## taxa table
timeStr,addedStr = populate_taxon_table(taxaList,session)
push_out(timeStr)
push_out(addedStr)

## gene table
timeStr,addedStr = populate_gene_table(taxaList,session)
push_out(timeStr)
push_out(addedStr)

## accession table
timeStr,addedStr = populate_accession_table(taxaList,session)
push_out(timeStr)
push_out(addedStr)

## go terms and go annotations tables
timeStr,addedStr = populate_go_tables(taxaList,session)
push_out(timeStr)
push_out(addedStr)

push_out("DATABASE - SUMMARY")
push_out("There are %s unique taxa "%session.query(Taxon).count())
push_out("There are %s unique genes   "%session.query(Gene).count())
push_out("There are %s unique accessions"%session.query(Accession).count())
print "\n"

fid.close()
