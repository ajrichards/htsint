import sys,os
import numpy as np
import matplotlib.pyplot as plt

from htsint import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation,db_connect

"""
Classes used to interact with gene ontology data
The data is stored in the htsint database

"""

class GeneOntology(object):
    "A class to interact with Gene Ontology data"
    
    def __init__(self,taxID,verbose=False):
        """
        Constructor
        
        taxID - is the NCBI taxid (int)
            i.e. for xenopus
            go = GeneOntology(8364)
        """

        ## variables
        self.taxID = taxID
        self.session,self.engine = db_connect(verbose=verbose)

        ## check if taxon is in database
        self.taxQuery = self.session.query(Taxon).filter_by(ncbi_id=self.taxID).first()
        if self.taxQuery == None:
            raise Exception("Taxon:%s not found in database"%self.taxID)

    def print_summary(self):
        """
        Print a summary of all Gene Ontology related information in database
        """
        print "\nGene Ontology Summary"
        print "----------------------------------"
        print "TaxID:       %s"%self.taxQuery.ncbi_id
        print "Species:     %s"%self.taxQuery.name
        print "Common Name: %s"%self.taxQuery.common_name_1
        
        geneQuery = self.session.query(Gene).filter_by(taxa_id=self.taxQuery.id)
        print "Num. Genes:  %s"%geneQuery.count()
        #print geneQuery[0]
        #print len(geneQuery)

        ## slow way
        for gene in geneQuery:
            annotQuery = self.session.query(GoAnnotation).filter_by(gene_id=gene.id).first()
            if annotQuery == None:
                continue
            print dir(annotQuery)
            #sys.exit()

        #print len(geneQuery)
        #print dir(self.taxQuery)
        print "\n"

if __name__ == "__main__":
    """
    Xenopus tropicalis - 8364
    Xenopus laevis     - 8355
    Mus musculus       - 10090
    """
    
    go = GeneOntology(8355)
    go.print_summary()
