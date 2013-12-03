import sys,os,cPickle
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
#from htsint import __basedir__
from htsint import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation,db_connect
from GeneOntologyLib import *

"""
Classes used to interact with gene ontology data
The data is stored in the htsint database
"""

class GeneOntology(object):
    "A class to interact with Gene Ontology data"
    
    def __init__(self,taxID=None,geneList=None,verbose=False):
        """
        Constructor
        
        taxID - is the NCBI taxid (int)
                i.e. for xenopus
                go = GeneOntology(8364)

                if no taxID is supplied then the class is in 'mixed' taxa mode
                and a geneList must be supplied.
        geneList - a list of ncbi gene ids

        """

        ## error checking
        if taxID == None and geneList == None:
            raise Exception("If taxid is not supplied a 'geneList' is required")

        ## start a database session
        self.session,self.engine = db_connect(verbose=verbose)

        if geneList != None:
            self.geneList = geneList
            self.taxID = None
        else:
            self.taxID = taxID
            self.check_taxon(taxID)
            self.taxQuery = self.session.query(Taxon).filter_by(ncbi_id=self.taxID).first()
            self.geneQuery = self.session.query(Gene).filter_by(taxa_id=self.taxQuery.id)
            self.geneList = [g.ncbi_id for g in self.geneQuery.all()]
            self.geneQuery = self.session.query(Gene).filter_by(taxa_id=self.taxQuery.id)
            self.annotQuery = self.session.query(GoAnnotation).filter_by(taxa_id=self.taxQuery.id)
        
    def check_taxon(self,taxID):
        """
        check if taxon is in database
        """

        taxQuery = self.session.query(Taxon).filter_by(ncbi_id=taxID).first()
        if taxQuery == None:
            raise Exception("Taxon:%s not found in database"%taxID)
        
    def print_summary(self):
        """
        Print a summary of all Gene Ontology related information in database
        """

        if self.taxID != None:
            print "\nSummary"
            print "----------------------------------"
            print "TaxID:       %s"%self.taxQuery.ncbi_id
            print "Species:     %s"%self.taxQuery.name
            print "Common Name: %s"%self.taxQuery.common_name_1
            print "Num. Genes:  %s"%self.geneQuery.count()
            print "Num. GO Annotations:  %s"%self.annotQuery.count()
            print "\n"
        else:
            print "\nSummary"
            print "----------------------------------"
            print "TaxID:       mixed taxa"
            print "Num. Genes:  %s"%len(self.geneList)
            print "\n"


    def get_terms(self,aspect='biological_process'):
        """
        get terms from the database associated with the gene list
        normally called from get_dicts
        """

        shortAspect = {"biological_process":"Process","molecular_function":"Function","cellular_component":"Component"}
        expEvidCodes = ["EXP","IDA","IPI","IMP","IGI","IEP"]
        compEvidCodes = ["ISS","ISO","ISA","ISM","IGC","RCA"]
        statEvidCodes = ["TAS","NAS","IC"]
        nonCuratedEvidCodes = ["IEA"]

        geneQueries = self.session.query(Gene).filter(Gene.ncbi_id.in_(self.geneList)).all()
        if len(geneQueries) != len(self.geneList):
            print "WARNING: %s/%s genes found in database"%(len(geneQueries),len(self.geneList))
        if len(geneQueries) == 0:
            raise Exception("There are no genes in gene list present in the database")

        gene2go = {}
        totalAnnotations = 0

        for gq in geneQueries:
            annotQuery = self.session.query(GoAnnotation).filter_by(gene_id=gq.id)
            
            if annotQuery.count() > 0:
                annotations = set([])
                for a in annotQuery.all():
                    termQuery = self.session.query(GoTerm).filter_by(id=a.go_term_id).first()
                    
                    ## filter by evidence code
                    if a.evidence_code in nonCuratedEvidCodes or a.evidence_code in compEvidCodes:
                        continue
                    
                    ## filter by aspect
                    if termQuery.aspect != shortAspect[aspect]:
                        continue

                    annotations.update([termQuery.go_id])
                if len(annotations) > 0:
                    gene2go[gq.ncbi_id] = list(annotations)
                    totalAnnotations += len(annotations)

        return gene2go

    def get_dicts(self,aspect='biological_process',filePath=None):
        """
        get the go2gene and gene2go dictionaries
        """

        if aspect not in ['biological_process','molecular_function','cellular_component']:
            raise Exception("Invalid aspect specified%s"%aspect)

        if filePath != None and os.path.exists(filePath):
            tmp = open(filePath,'r')
            gene2go,go2gene = cPickle.load(tmp)
            tmp.close()
            return gene2go,go2gene

        ## gene2go
        print "...creating gene2go dictionary -- this may take several minutes or hours depending on the number of genes"
        gene2go = self.get_terms()

        ## go2gene
        print "...creating go2gene dictionary -- this may take several minutes"
        go2gene = {}
        for gene,terms in gene2go.iteritems():
            for term in terms:
                if go2gene.has_key(term) == False:
                    go2gene[term] = set([])
                go2gene[term].update([gene])

        for term,genes in go2gene.iteritems():
            go2gene[term] = list(genes)
        
        if filePath == None:
            print "INFO: For large gene list it is useful to specify a file path for the pickle file"
        else:
            tmp = open(filePath,'w')
            cPickle.dump([gene2go,go2gene],tmp)
            tmp.close()

        return gene2go,go2gene

    def create_gograph(self,aspect='biological_process',termsPath=None,graphPath=None):
        """
        creates the go graph for a given species
        aspect = 'biological_process','molecular_function' or 'cellular_component'

        A go graph can be created with any gene list
        The genes have to be present in the database
        """

        ## error checking
        if aspect not in ['biological_process','molecular_function','cellular_component']:
            raise Exception("Invalid aspect specified%s"%aspect)

        ## load pickled version if present
        if graphPath != None and os.path.exists(graphPath):
            G = nx.read_gpickle(graphPath)
            return G

        ## load all the ontology related information as a set of dictionaries
        _goDict = read_ontology_file()    
        goDict = _goDict[aspect]
        gene2go,go2gene = self.get_dicts(filePath=termsPath)

        print "...creating go term graph -- this may take several minutes or hours depending on the number of genes"
        
        ## calculate the information content for each term -ln(p(term))
        edgeDict = self.get_weights_by_ic(goDict,go2gene)

        ## initialize the graph
        G = nx.DiGraph()
        for nodes,weight in edgeDict.iteritems():
            parent,child = nodes.split("#")
            G.add_edge(parent,child,weight=weight)
        
        ## save it to pickle format
        if graphPath != None:
            nx.write_gpickle(G, graphPath)
            print '...saving pickle graph'

        return G

    def get_weights_by_ic(self,goDict,go2gene):
        """
        get the weight of a term using information content
        returns a networkx edge dictionary
        """

        print "...... terms", len(go2gene.keys())

        total = 0
        for term,genes in go2gene.iteritems():
            total += len(genes)

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
        return edgeDict


if __name__ == "__main__":
    """
    Xenopus tropicalis      - 8364
    Xenopus laevis          - 8355
    Mus musculus            - 10090
    Drosophila melanogaster - 7227
    """

    ## testing with a species
    go = GeneOntology(taxID=7227)
    go.print_summary()
    #go.create_gograph()

    ## testing with a gene list
    geneList = ['30970','30971','30972','30973','30975','30976','30977','30978','30979','30980',
                '30981','30982','30983','30984','30985','30986','30988','30990','30991','30994',
                '30995','30996','30998','31000','31001','31002','31003','31004','31005','31006']

    go = GeneOntology(geneList=geneList)
    go.print_summary()
    go.create_gograph()

    #go.create_gograph()
