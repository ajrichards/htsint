#!/usr/bin/env python
"""
A class that prints and generates summary information about a species list
with respect to the information contained in the database
"""

__author__ = "Adam Richards"

import csv,re,sys
from htsint.database import db_connect,Taxon,Gene,GoAnnotation

class TaxaSummary(object):
    """
    Summarize the information in a taxa list
    """

    def __init__(self,taxaFilePath,delimiter="\t",upass=''):
        """
        The input file is a tab delimited file 
        The delimiter can also be specified
        The first column is the ncbi taxa id
        All additional columns are optional and user specified 
        """

        self.session,self.engine = db_connect(upass=upass)
        self.taxaList = self.read_file(taxaFilePath)
        print self.taxaList

        ## get the taxa ids associated wtih the taxa list
        taxaIds = [tid.id for tid in self.session.query(Taxon).filter(Taxon.ncbi_id.in_(self.taxaList)).all()]   
        print taxaIds

        blah = [tid for tid in self.session.query(Taxon).filter(Taxon.ncbi_id.in_(self.taxaList)).all()]
        for b in blah:
            print dir(b)
            print b.name,b.ncbi_id,b.common_name_1

        sys.exit()


        annots = {}

        for t,taxaId in enumerate(taxaIds):
            print t,taxaId
            annots[self.taxaList[t]] = {'bytaxa':self.session.query(GoAnnotation).filter_by(taxa_id=taxaId).all()}




    def read_file(self,taxaFilePath):
        """
        Use csv to read in the species 
        """
        
        fid = open(taxaFilePath,'rU')
        reader = csv.reader(fid,delimiter='\t')
        taxaList = []

        for linja in reader:
            if re.search("\D",linja[0]):
                continue
            taxaList.append(linja[0])

        fid.close()

        return taxaList


    



if __name__ == "__main__":
    print "Running..."
