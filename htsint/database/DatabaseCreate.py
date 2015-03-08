#!/usr/bin/env python
"""
creates htsint database
(1) Edit/create the configure.py file and specify the database parameters
(2) Run 'DatabaseFetch.py'
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
from htsint import Configure
from DatabaseTables import Base,Taxon,Gene,Uniprot,GoTerm,GoAnnotation
from DatabaseTools import db_connect, get_file_sizes,print_db_summary
from DatabaseTools import populate_taxon_table,populate_gene_table,populate_uniprot_table
from DatabaseTools import populate_go_terms, populate_go_annotations
from GeneOntologyLib import read_annotation_file,get_annotation_file,get_total_annotations

class DatabaseCreate(object):
    """
    Database population class
    """

    def __init__(self):
        """
        Constructor
        """

        self.config = Configure()

        ## ensure config is setup 
        for key in ['data','dbname']:
            if self.config.log[key] == '':
                raise Exception("You must modify the config file before running DatabaseFetch.py")

        dataDir = self.config.log['data']

        if not os.path.isdir(dataDir):
            raise Exception("Specified htsint data directory does not exist %s"%dataDir)

        self.taxaList = self.config.log['taxa']

    def run(self):
        ## prepare a log file
        fid = open(os.path.join(self.config.log['data'],'createdb.log'),'wa')
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

        ## determine file sizes
        timeStart = time.time()
        push_out("determining filesizes...")
        idmapCount,geneInfoCount = get_file_sizes()

        push_out("extracting taxa list...")
        totalAnnotations = get_total_annotations()
        push_out("...extraction time: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

        ## taxa table
        push_out("Populating the database taxa table")
        timeStr,addedStr = populate_taxon_table(engine)
        push_out(timeStr)
        push_out(addedStr)

        ## gene table
        push_out("Populating the database with %s genes"%(geneInfoCount))
        timeStr,addedStr = populate_gene_table(geneInfoCount,session,engine)
        push_out(timeStr)
        push_out(addedStr)
        
        ## uniprot table
        push_out("Populating the database with %s uniprot entries"%(idmapCount))
        timeStr,addedStr = populate_uniprot_table(idmapCount,session,engine)
        push_out(timeStr)
        push_out(addedStr)

        ## populate the go-terms
        push_out("Populating the database with for go terms...")
        timeStr,addedStr = populate_go_terms(engine)
        push_out(timeStr)
        push_out(addedStr)

        ## populate the go-annotations
        push_out("Populating the database with for go annotations...")
        timeStr,addedStr,ignored = populate_go_annotations(totalAnnotations,session,engine)
        push_out(timeStr)
        push_out(addedStr)
        push_out("There were %s uniprot annotations ignored"%str(ignored[0]))
        push_out("There were %s gene annotations ignored"%str(ignored[1]))

        print_db_summary()
        fid.close()

if __name__ == "__main__":
    print "Running..."
    dbc = DatabaseCreate()
    dbc.run()
