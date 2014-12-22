import sys,os,cPickle
import numpy as np
import networkx as nx
from sqlalchemy.sql import select
from htsint.database import Base,Taxon,Gene,Uniprot,GoTerm,GoAnnotation,db_connect
from htsint.database import read_ontology_file,fetch_taxa_annotations
from htsint.stats import EmpiricalCdf

"""
Classes used to interact with gene ontology data
The data is stored in the htsint database
"""

class GeneOntology(object):
    "A class to interact with Gene Ontology data"
    
    def __init__(self,taxaList,verbose=False,upass='',idType='ncbi',useIea=True,\
                 aspect='biological_process'):
        """
        Constructor
        
        taxaList a list of NCBI taxa ids
        """

        ## error checking
        idType = idType.lower()
        if idType not in ['uniprot','ncbi']:
            raise Exception("Invalid idType argument in fetch annotations use 'uniprot' or 'ncbi'")

        ## start a database session
        self.session,self.engine = db_connect(verbose=verbose,upass=upass)

        ## global variables
        self.taxaList = taxaList
        self.idType = idType
        self.useIea = useIea
        self.aspect = aspect

    def summarize(self,refTaxon,termsPath):
        """
        GO object summary and sanity check
        """

        refTaxon = str(refTaxon)
        if refTaxon not in self.taxaList:
            raise Exception("refTaxon not present in taxaList")

        conn = self.engine.connect()
        gene2go,go2gene = self.get_dicts(termsPath=termsPath)

        s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(self.taxaList))
        _taxaQueries = conn.execute(s)
        taxaQueries = _taxaQueries.fetchall()
        taxaMap = dict([(str(r['ncbi_id']),str(r['id'])) for r in taxaQueries])

        gene2id = {}
        for tquery in taxaQueries: 
            s = select([Gene.id,Gene.ncbi_id],Gene.taxa_id==tquery['id'])
            _geneQueries = conn.execute(s)
            taxaDict = dict([(str(r['ncbi_id']),str(r['id'])) for r in _geneQueries.fetchall()])
            if str(tquery['ncbi_id']) == refTaxon:
                refGenes = taxaDict.copy()

            print("there are  %s genes from %s (%s)"%(len(taxaDict.keys()),tquery['name'],tquery['ncbi_id']))
            gene2id.update(taxaDict)

        ## check for unmatched genes
        unmatched = 0
        for gene in gene2go.iterkeys():
            if not gene2id.has_key(gene):
                unmatched += 1
        
        print("Summary")
        print("IEA annotations: %s"%self.useIea)
        print("total genes in combined taxa: %s"%(len(gene2id.keys())))
        if unmatched > 0:
            print("WARNING: there were unmatched genes unmatched: %s"%unmatched)
        print("total genes with at least one annotation: %s"%(len(gene2go.keys())))
        print("total unique annotations: %s"%(len(go2gene.keys())))
        print("---------------------")

        _gene2go,_prot2go = fetch_taxa_annotations([refTaxon],self.engine,aspect=self.aspect,\
                                                   useIea=self.useIea)

        total= 0
        for k,v in _gene2go.iteritems():
            total += len(v)
        print('RefTaxa genes: %s'%(len(refGenes.keys())))
        print('Only RefTaxa: %s annotated genes, %s total annotations'%(len(_gene2go.keys()), total))
        
        total= 0
        for k,v in gene2go.iteritems():
            total += len(v)
        print('With additional taxa: %s annotated genes, %s total annotations'%(len(gene2go.keys()), total))
        print('Percent annotation: %s'%(float(len(gene2go.keys())) / float(len(refGenes.keys()))))
                                       

    def check_taxon(self,taxID):
        """
        check if taxon is in database
        """

        taxQuery = self.session.query(Taxon).filter_by(ncbi_id=taxID).first()
        if taxQuery == None:
            raise Exception("Taxon:%s not found in database"%taxID)
        
    def get_dicts(self,termsPath=None,log=None):
        """
        get the go2gene and gene2go dictionaries
        log - csv.writer object
        """

        conn = self.engine.connect()

        if self.aspect not in ['biological_process','molecular_function','cellular_component']:
            raise Exception("Invalid aspect specified%s"%self.aspect)

        if termsPath != None and os.path.exists(termsPath):
            tmp = open(termsPath,'r')
            gene2go,go2gene = cPickle.load(tmp)
            tmp.close()
            return gene2go,go2gene

        ## gene2go
        print "...creating gene2go dictionary -- this may take several minutes or longer on the number of genes"
        gene2go,prot2go = fetch_taxa_annotations(self.taxaList,self.engine,aspect=self.aspect,\
                                                 useIea=self.useIea)

        print "...creating go2gene dictionary -- this may take several minutes"
        go2gene = {}
        for gene,terms in gene2go.iteritems():
            for term in terms:
                if go2gene.has_key(term) == False:
                    go2gene[term] = set([])
                go2gene[term].update([gene])

        for term,genes in go2gene.iteritems():
            go2gene[term] = list(genes)
        
        if termsPath == None:
            print "INFO: For large gene list it is useful to specify a file path for the pickle file"
        else:
            tmp = open(termsPath,'w')
            cPickle.dump([gene2go,go2gene],tmp)
            tmp.close()

        return gene2go,go2gene

    def create_gograph(self,termsPath=None,graphPath=None):
        """
        creates the go graph for a given species
        aspect = 'biological_process','molecular_function' or 'cellular_component'

        A go graph can be created with any gene list
        The genes have to be present in the database
        """

        print '...creating gograph'
        ## error checking
        if self.aspect not in ['biological_process','molecular_function','cellular_component']:
            raise Exception("Invalid aspect specified%s"%self.aspect)

        ## load pickled version if present
        if graphPath != None and os.path.exists(graphPath):
            G = nx.read_gpickle(graphPath)
            return G

        ## load all the ontology related information as a set of dictionaries
        _goDict = read_ontology_file()
        goDict = _goDict[self.aspect]
        gene2go,go2gene = self.get_dicts(termsPath=termsPath)

        #for a,b in gene2go.iteritems():
        #    print a,b

        print "...creating go term graph -- this may take several minutes or hours depending on the number of genes"
        ## calculate the term-term edges using term IC (-ln(p(term)))
        edgeDict = self.get_weights_by_ic(goDict,go2gene)

        if len(edgeDict) == 0:
            raise Exception("No edges were found -- annotation information may be two sparse")

        ## terms without distances set to max
        print "...finding term-term distances"
        allDistances = np.array(edgeDict.values())
        maxDistance = allDistances.max()
        for parent,children in goDict.iteritems():
            for child in list(children):
                if edgeDict.has_key(parent+"#"+child) or edgeDict.has_key(child+"#"+parent):
                    continue
                edgeDict[parent+"#"+child] = maxDistance

        ## add the edges through shared genes
        print '...adding term-term edges through shared genes'
        eCDF = EmpiricalCdf(allDistances)
        p5 = eCDF.get_percentile(0.05) 
        newEdges = 0
        for termI,genesI in go2gene.iteritems():
            for termJ,genesJ in go2gene.iteritems():

                if edgeDict.has_key(termI+"#"+termJ) or edgeDict.has_key(termJ+"#"+termI):
                    continue

                sharedGenes = len(list(set(genesI).intersection(set(genesJ))))
                if sharedGenes == 0:
                    continue

                ## add new edge in both directions
                edgeDict[termI+"#"+termJ] = float(p5) / float(sharedGenes)
                newEdges += 1

        ## initialize the graph
        G = nx.Graph()
        for nodes,weight in edgeDict.iteritems():
            parent,child = nodes.split("#")
            G.add_edge(parent,child,weight=weight)
        
        ## save mst to pickle format
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
        
        return edgeDict

if __name__ == "__main__":
    """
    Xenopus tropicalis      - 8364
    Xenopus laevis          - 8355
    Mus musculus            - 10090
    Drosophila melanogaster - 7227
    """

    ## testing with a species
    go = GeneOntology(['7091'],useIea=True)
    go.create_gograph()

