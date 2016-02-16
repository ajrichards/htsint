import sys,os,cPickle
import numpy as np
import networkx as nx
from sqlalchemy.sql import select
from htsint.database import Base,Taxon,Gene,Uniprot,GoTerm,GoAnnotation,db_connect
from htsint.database import read_ontology_file,fetch_taxa_annotations

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
        gene2go,go2gene = self.load_dicts(termsPath=termsPath)

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
        

    def create_dicts(self,termsPath,accepted=None):
        """
        get the go2gene and gene2go dictionaries
        'accepted' - list of genes that restrict included terms to a particular list
        """

        conn = self.engine.connect()

        ## error checking
        if self.aspect not in ['biological_process','molecular_function','cellular_component']:
            raise Exception("Invalid aspect specified%s"%self.aspect)

        ## gene2go
        print("...creating gene2go dictionary -- this may take several minutes or longer depending on the number of genes")
        _gene2go,prot2go = fetch_taxa_annotations(self.taxaList,self.engine,aspect=self.aspect,\
                                                 useIea=self.useIea)

        print("...creating go2gene dictionary -- this may take several minutes")
        go2gene = {}
        gene2go = {}
        for gene,terms in _gene2go.iteritems():
            if accepted and gene not in accepted:
                continue

            gene2go[gene] = terms
            for term in terms:
                if go2gene.has_key(term) == False:
                    go2gene[term] = set([])
                go2gene[term].update([gene])

        for term,genes in go2gene.iteritems():
            go2gene[term] = list(genes)
        
        ## pickle the dictionaries    
        tmp = open(termsPath,'w')
        cPickle.dump([gene2go,go2gene],tmp)
        tmp.close()

    def load_dicts(self,termsPath=None,log=None):
        """
        get the go2gene and gene2go dictionaries
        """

        if self.aspect not in ['biological_process','molecular_function','cellular_component']:
            raise Exception("Invalid aspect specified%s"%self.aspect)

        if termsPath != None and os.path.exists(termsPath):
            tmp = open(termsPath,'r')
            gene2go,go2gene = cPickle.load(tmp)
            tmp.close()
            return gene2go,go2gene

        return None,None

    def create_gograph(self,termsPath=None,graphPath=None,addShared=False):
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
        gene2go,go2gene = self.load_dicts(termsPath=termsPath)
        if gene2go == None or go2gene == None:
            raise Exception("Could not load gene2go or go2Gene-- do they exist?")

        if len(go2gene.keys()) < 5:
            raise Exception("There are insufficient terms present to create a graph %s"%(len(go2gene.keys())))

        _goDict = read_ontology_file()
        goDict = _goDict[self.aspect]
        
        print "...creating go term graph -- this may take several minutes or hours depending on the number of genes"
        ## calculate the term-term edges using term IC (-ln(p(term)))
        edgeDict = self.get_weights_by_ic(goDict,go2gene)

        if len(edgeDict) == 0:
            raise Exception("No edges were found -- annotation information may be two sparse")

        ## terms without distances set to max
        print("...finding term-term distances")
        allDistances = np.array(edgeDict.values())
        maxDistance = allDistances.max()
        for parent,children in goDict.iteritems():
            for child in list(children):
                if edgeDict.has_key(parent+"#"+child) or edgeDict.has_key(child+"#"+parent):
                    continue
                edgeDict[parent+"#"+child] = maxDistance

        ## add the edges through shared genes
        if addShared == True:
            print('...adding term-term edges through shared genes')
            p5 = np.percentile(allDistances,0.05)
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
        
        def trim_leaf_nodes(G):
            """
            trim nodes with degree = 1 that are not present in go2gene
            """

            degreeDict = nx.degree(G)
            removed = 0
            for node,degree in degreeDict.iteritems():

                if degree > 1:
                    continue

                if go2gene.has_key(node):
                    continue

                removed += 1
                G.remove_node(node)

            return G,removed

        ## trim unecessary leaf nodes
        print("trimming leaf nodes...")
        removed = 1
        print("nodes, edges")
        print("%s, %s"%(G.number_of_nodes(), G.number_of_edges()))
        while removed > 0:
            G,removed = trim_leaf_nodes(G)
            print("%s, %s"%(G.number_of_nodes(), G.number_of_edges()))

        ## trim nodes with degree 2
        print("trimming degree=2 nodes...")
        toRemove = []
        toAdd    = []
        degreeDict = nx.degree(G)
        for node,degree in degreeDict.iteritems():

            if degree != 2:
                continue

            if go2gene.has_key(node):
                continue

            neighbors = [x for x in nx.all_neighbors(G,node)]

            if neighbors[0] in toRemove or neighbors[1] in toRemove:
                continue

            edge1 = G.get_edge_data(neighbors[0],node)['weight']
            edge2 = G.get_edge_data(neighbors[1],node)['weight']

            toRemove.append(node)
            toAdd.append((neighbors[0],neighbors[1],np.array([edge1,edge2]).mean()))

        for node in toRemove:
            G.remove_node(node)
        for edge in toAdd:
            G.add_edge(edge[0],edge[1],weight=edge[2])

        ## trim unecessary leaf nodes
        print("trimming leaf nodes...")
        removed = 1
        print("nodes, edges")
        print("%s, %s"%(G.number_of_nodes(), G.number_of_edges()))
        while removed > 0:
            G,removed = trim_leaf_nodes(G)
            print("%s, %s"%(G.number_of_nodes(), G.number_of_edges()))

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

