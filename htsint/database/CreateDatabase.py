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
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation
from DatabaseTools import db_connect
from DatabaseTools import populate_taxon_table,populate_gene_table,populate_accession_table,populate_go_tables

from GeneOntologyLib import read_idmapping_file

## prepare a log file
fid = open('createdb.log','w')
writer = csv.writer(fid)

def push_out(line):
    writer.writerow([line])
    print line

push_out(sys.argv[0])
push_out(time.asctime())
push_out("Getting ready to create database...")

## conect to the database
session,engine = db_connect(verbose=False)

## create the tables (uncomment to erase everything first)
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine) 

## read the annotation file
taxaList,annotations = read_idmapping_file()

## taxa table
print "total", len(taxaList)
taxaList = taxaList[:1]
print taxaList
push_out("Attempting to populate the database with %s taxa"%(len(taxaList)))
timeStr,addedStr = populate_taxon_table(taxaList,session)
push_out(timeStr)
push_out(addedStr)

## gene table
timeStr,addedStr = populate_gene_table(taxaList,annotations,session)
push_out(timeStr)
push_out(addedStr)


print 'here we go'

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
