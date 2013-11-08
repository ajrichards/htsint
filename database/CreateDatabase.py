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
    Use PrimeDatabase.py to further populate the database

For more information about interacting with databases using SQLalchemy
    http://docs.sqlalchemy.org/en/rel_0_8/orm/tutorial.html


Because gene names are added/removed/changed on a regular basis
It is easier when updating to re-run this script and then follow up
with the PopulateDatabase.py afterwards. 

"""

### make imports
import sys,os,re,time,csv
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation
from DatabaseTools import db_connect
from DatabaseTools import populate_taxon_table,populate_gene_table,populate_accession_table,populate_go_tables
from config import CONFIG

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

## taxon table
taxaList = ["7227"] + CONFIG['taxa']

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
