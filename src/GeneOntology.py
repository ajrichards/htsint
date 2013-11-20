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

    def create_gograph(self,aspect='biological_process',geneList=None):
        """
        creates the go graph for a given species
        aspect = 'biological_process','molecular_function' or 'cellular_component'
        geneList = is by default found using the database, however it may also be specified

        A go graph can be created with any gene list
        The genes have to be present in the database
        """

        shortAspect = {"biological_process":"Process","molecular_function":"Function","cellular_component":"Component"}
        expEvidCodes = ["EXP","IDA","IPI","IMP","IGI","IEP"]
        compEvidCodes = ["ISS","ISO","ISA","ISM","IGC","RCA"]
        statEvidCodes = ["TAS","NAS","IC"]
        nonCuratedEvidCodes = ["IEA"]
        _goDict = read_ontology_file()    
        goDict = _goDict['biological_process']

        ## error checking
        if aspect not in ['biological_process','molecular_function','cellular_component']:
            raise Exception("Invalid aspect specified%s"%aspect)

        ## get the gene list
        if geneList == None:
            geneList = [g.symbol for g in self.geneQuery.all()]

        ## get the terms associated with each gene
        go2gene = {}
        for item in self.annotQuery.all():
            ## include only certain evidence codes
            if item.evidence_code  in nonCuratedEvidCodes or item.evidence_code in compEvidCodes:
                continue

            ## include only the correct aspect
            termQuery = self.session.query(GoTerm).filter_by(id=item.go_term_id).first()
            if termQuery.aspect != shortAspect[aspect]:
                continue
            
            ## update the dictionary
            symbol = self.geneQuery.filter_by(id=item.gene_id).first().symbol
            if go2gene.has_key(termQuery.go_id) == False:
                go2gene[termQuery.go_id] = set([])
            go2gene[termQuery.go_id].update([symbol])
        
        ## calculate the information content for each term -ln(p(term))
        total = 0
        for term,genes in go2gene.iteritems():
            total += len(list(genes))

        edgeDict = {}
        for parent,children in goDict.iteritems():
            if go2gene.has_key(parent) == False:
                continue

            parentIC = -np.log(float(len(list(go2gene[parent])))  / float(total))
            for child in list(children):
                if go2gene.has_key(child) == False:
                    continue

                childIC = -np.log(float(len(list(go2gene[child])))  / float(total))
                distance = np.abs(parentIC - childIC)
                edgeDict[parent+"#"+child] = distance

        print 'total annotations', total
        print 'done'

        ## initialize the graph
        G = nx.DiGraph()
        for nodes,weight in edgeDict.iteritems():
            parent,child = nodes.split("#")
            G.add_edge(parent,child,weight=weight)
        
        ## save it to pickle format
        filePath = os.path.join(__basedir__,'graphs','gograph_%s.pickle'%(self.taxQuery.ncbi_id))
        nx.write_gpickle(G, filePath)
        print 'saving pickle graph'                      


    def get_weight_by_ic(term):
        """
        get the weight of a term
        """


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
