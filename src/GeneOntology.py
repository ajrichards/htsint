import sys,os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from htsint import __basedir__
from htsint import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation,db_connect
from GeneOntologyLib import *

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

        ## other database queries
        self.geneQuery = self.session.query(Gene).filter_by(taxa_id=self.taxQuery.id)
        self.annotQuery = self.session.query(GoAnnotation).filter_by(taxa_id=self.taxQuery.id)
        
    def print_summary(self):
        """
        Print a summary of all Gene Ontology related information in database
        """
        print "\nGene Ontology Summary"
        print "----------------------------------"
        print "TaxID:       %s"%self.taxQuery.ncbi_id
        print "Species:     %s"%self.taxQuery.name
        print "Common Name: %s"%self.taxQuery.common_name_1
        print "Num. Genes:  %s"%self.geneQuery.count()
        print "Num. GO Annotations:  %s"%self.annotQuery.count()
        print "\n"

    def create_gograph(self):


        read_ontology_file()

        ## check to see if ontology file is present
        


        #edgeDict = {"w4a":("X1","X4"),"w5a":("X1","X5"),"w6a":("X1","X6"),"w7a":("X1","X7"),
        #            "w4b":("X2","X4"),"w5b":("X2","X5"),"w6b":("X2","X6"),"w7b":("X2","X7"),
        #            "w4c":("X3","X4"),"w5c":("X3","X5"),"w6c":("X3","X6"),"w7c":("X3","X7"),
        #            "w1":("X1","T"),"w2":("X2","T"),"w3":("X3","T"),"w4":("X4","T"),
        #            "w5":("X5","T"),"w6":("X6","T"),"w7":("X7","T")}

        ## initialize the graph
        #G = nx.Graph()
        #for node in ["X1","X2","X3","X4","X5","X6","X7","T"]:
        #    G.add_node(node)
        #
        #for edgeName,edge in edgeDict.iteritems():
        #    G.add_edge(edge[0],edge[1],weight=edgeWeights[edgeName])



if __name__ == "__main__":
    """
    Xenopus tropicalis      - 8364
    Xenopus laevis          - 8355
    Mus musculus            - 10090
    Drosophila melanogaster - 7227
    """
    
    go = GeneOntology(7227)
    go.print_summary()
    go.create_gograph()
