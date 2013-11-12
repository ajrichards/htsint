#!/usr/bin/env python
"""
See CreateDatabase.py before running DatabaseAppend.py

DatabaseAppend.py takes as input one or more taxon ids
"""

### make imports
import sys,os,re,time,csv,getopt
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation
from DatabaseTools import db_connect
from DatabaseTools import populate_taxon_table,populate_gene_table,populate_accession_table,populate_go_tables

## read in input file
if len(sys.argv) < 1:
    print sys.argv[0] + " -t taxa_list"
    sys.exit()

try:
    optlist, args = getopt.getopt(sys.argv[1:], 't:')
except getopt.GetoptError:
    print sys.argv[0] + " -t taxa_list"
    sys.exit()

taxaList = None
for o, a in optlist:
    if o == '-t':
        taxaList = a

## error checking
if taxaList == None:
    print "\nINPUT ERROR: incorrect arguments"
    print sys.argv[0] + " -t taxa_list"
    print "For example..."
    print sys.argv[0] + " -t 8355,8364\n"
    sys.exit()

print taxaList
_taxaList = taxaList.split(",")

## prepare a log file
print("Connecting to the database...")

## conect to the database
session,engine = db_connect(verbose=False)
print("Appending to database...")

taxaList = []
for taxID in _taxaList:
    query = session.query(Taxon).filter_by(ncbi_id=taxID).first()
    if query == None:
        taxaList.append(taxID)
    else:
        print("The taxon %s is already present in the database skipping..."%taxID)

if len(taxaList) == 0:
    print "All taxa already present in database"
    sys.exit()

print("adding... %s"%str(taxaList))

## taxon table
timeStr,addedStr = populate_taxon_table(taxaList,session)

## gene table
timeStr,addedStr = populate_gene_table(taxaList,session)

## accession table
timeStr,addedStr = populate_accession_table(taxaList,session)

## go terms and go annotations tables
timeStr,addedStr = populate_go_tables(taxaList,session)


print("DATABASE - SUMMARY")
print("There are %s unique taxa "%session.query(Taxon).count())
print("There are %s unique genes   "%session.query(Gene).count())
print("There are %s unique accessions"%session.query(Accession).count())
print "\n"
