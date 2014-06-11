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
        self.taxaList = list(set(self.read_file(taxaFilePath)))
        print self.taxaList

        ## get the taxa ids associated wtih the taxa list
        taxaRows = []
        for taxonId in self.taxaList:
            query = self.session.query(Taxon).filter_by(ncbi_id=taxonId).all()
            if len(query) > 1:
                raise Exception("Taxon id '%s' return multiple values from database"%taxonId)
            if query == None:
                raise Exception("Invalid taxon '%s' provided in list"%taxonId)
            taxaRows.append(query[0])

        ## print taxa name information
        print("------------------")
        print("taxa names")
        print("------------------")
        for tr in taxaRows:
            print("id=%s,name=%s,common_name=%s"%(tr.ncbi_id,tr.name,tr.common_name_1))
    
        ## print gene infomation
        print("------------------")
        print("gene info")
        print("------------------")
        for taxaRow in taxaRows:
            queryGenes = self.session.query(Gene).filter_by(taxa_id=taxaRow.id).all()
            print("id=%s,num_genes=%s"%(tr.ncbi_id,len(queryGenes)))
            



        
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
