#!/usr/bin/env python
"""
Used to test the database
"""

### make imports
import sys,os,re

from config import CONFIG
from DatabaseTables import Base,Taxon,Gene,Accession
from DatabaseTools import db_connect

try:
    from sqlalchemy_schemadisplay import create_schema_graph
    createGraph = True
except:
    createGraph = False

## conect to the database                                                                                                                                                         
session,engine = db_connect(verbose=False)

## test the 'Taxon' table
query = session.query(Taxon).filter_by(ncbi_id='7227').first() 

if query == None:
    print "ERROR: init taxon not found."
    sys.exit()

if query.ncbi_id != 7227:
    print "ERROR: Bad match to taxon id"
    print query.ncbi_id

if query.name != "Drosophila melanogaster":
    print "ERROR: Bad match to taxon name"
    print query.name
    sys.exit()

if "fruit fly" not in [query.common_name_1, query.common_name_2, query.common_name_3]:
    print "ERROR: Bad match to common name"
    print query.name
    sys.exit()

## test the 'Gene' table
if session.query(Gene).count() < 4:
    print "ERROR: not enough genes found"
    sys.exit()

query = session.query(Gene).filter_by(ncbi_id='3771877').first() 

if query.symbol != 'Adh':
    print "ERROR: bad gene symbol found"
    print query.symbol
    sys.exit()

## test the 'Accession' table
if session.query(Accession).count() < 100:
    print "ERROR: not enough accessions found"
    sys.exit()

if createGraph == True:
    # create the pydot graph object by autoloading all tables via a bound metadata object                                                                                             
    graph = create_schema_graph(metadata=Base.metadata,
                                show_datatypes=False,   # can get large with datatypes                                                                                                
                                show_indexes=False,     # ditto for indexes                                                                                                           
                                rankdir='LR',           # From left to right (instead of top to bottom)                                                                               
                                concentrate=False       # Don't try to join the relation lines together                                                                               
                            )
    #graph.write_svg('dbschema.svg')                     # write out the file 
    graph.write_png('dbschema.png')                     # write out the file 
else:
    print "Not creating schema figure because 'sqlalchemy_schemadisplay' is not installed"

print 'all tests pass.'
