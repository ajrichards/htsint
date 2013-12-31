#!/usr/bin/env python
"""
See CreateDatabase.py before running DatabaseAppend.py
DatabaseAppend takes as input one or more taxon ids

it may be run as a script:

    ~$ python DatabaseAppend -t 8355,8364

or from withing python

    >>> from htsint import DatabaseAppend
    >>> da = DatabasesAppend([8355,8364])
    >>> da.run()

"""

### make imports
import sys,os,re,time,csv,getopt
from DatabaseTools import db_connect
from DatabaseTools import Taxon,populate_taxon_table,populate_gene_table,populate_accession_table,populate_go_tables

class DatabaseAppend(object):
    """
    appends taxa to the existing database
    """
    
    def __init__(self,taxa,resultsDir="."):

        ## error checking
        if type(taxa) != type([]) or len(taxa) == 0:
            raise Exception("Invalid taxa list provided - must be of type list and must have ids")

        ## conect to the database
        self.session,self.engine = db_connect(verbose=False)
        print("Appending to database...")

        self.taxaList = list(set(taxa))

    def is_valid_list(self):
        """
        checks the taxa list before running
        """

        taxaList = []
        for taxID in self.taxaList:
            if not re.search("\d",str(taxID)) or re.search("\D",str(taxID)):
                print("The taxon %s is not a valid one skipping..."%taxID)
                continue
            query = self.session.query(Taxon).filter_by(ncbi_id=taxID).first()
            if query == None:
                taxaList.append(taxID)
            else:
                print("The taxon %s is already present in the database skipping..."%taxID)

        print 'list', taxaList
        self.taxaList = taxaList

        if len(self.taxaList) == 0:
            print "All taxa already present in database"
            return False
        else:
            return True
    
    def run(self):
        """
        runs DatabaseAppend
        """
        
        goFlag = self.is_valid_list()
        if goFlag == False:
            return

        print("adding... %s new ones"%len(self.taxaList))

        ## taxon table
        timeStr,addedStr = populate_taxon_table(self.taxaList,self.session)

        ## gene table 
        timeStr,addedStr = populate_gene_table(self.taxaList,self.session)

        ## accession table
        timeStr,addedStr = populate_accession_table(self.taxaList,self.session)

        ## go terms and go annotations tables
        timeStr,addedStr = populate_go_tables(self.taxaList,self.session)

        print "database append complete."


if __name__ == "__main__":
    
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

    taxa = taxaList.split(",")
    da = DatabaseAppend(taxa)
    da.run()
