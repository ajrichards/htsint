#!/usr/bin/env python
"""
Used to test the database

The schema display is borrowed from:
http://www.sqlalchemy.org/trac/wiki/UsageRecipes/SchemaDisplay
"""

### make imports
import sys,os,re
from DatabaseTables import Base,Taxon,Gene,Uniprot
from DatabaseTools import db_connect

try:
    from sqlalchemy_schemadisplay import create_schema_graph
    createGraph = True
except:
    createGraph = False

## conect to the database                                                            
session,engine = db_connect(verbose=False)

## test the 'Taxon' table
testID = '7227'
query = session.query(Taxon).filter_by(ncbi_id=testID).first() 

if query == None:
    print "ERROR: init taxon not found."
    sys.exit()

if int(query.ncbi_id) != int(testID):
    print "ERROR: Bad match to taxon id"
    print query.ncbi_id,testID

if query.name != "Drosophila melanogaster":
    print "ERROR: Bad match to taxon name"
    print query.name
    sys.exit()

if "fruit fly" not in [query.common_name_1, query.common_name_2, query.common_name_3]:
    print "ERROR: Bad match to common name"
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

## test the 'Uniprot' table
if session.query(Uniprot).count() < 100:
    print "ERROR: not enough uniprot found"
    sys.exit()

if createGraph == True:
    # create the pydot graph object by autoloading all tables via a bound metadata object
    graph = create_schema_graph(metadata=Base.metadata,
                                show_datatypes=False,   # can get large with datatypes
                                show_indexes=False,     # ditto for indexes
                                rankdir='LR',           # From left to right (instead of top to bottom)
                                concentrate=False       # Don't try to join the relation lines together 
                            )
    graph.write_svg('dbschema.svg')                     # write out the file 
    #graph.write_png('dbschema.png')                     # write out the file 
else:
    print "Not creating schema figure because 'sqlalchemy_schemadisplay' is not installed"

print 'all tests pass.'
