#!/usr/bin/env python
"""
A class to represent and draw gene sets
"""

__author__ = "Adam Richards"

import os,cPickle,csv
import networkx as nx
import matplotlib.pyplot as plt
from htsint.database import db_connect,Gene,Taxon,GoTerm
from htsint.blast import BlastMapper

class GeneSet(object):
    """
    gene set class
    """

    def __init__(self,genesetFile,gene2go):
        """
        Constructor
        geneList - list
        gene2go - dictionary

        """

        ## setup db
        self.session, self.engine = db_connect()
        self.conn = self.engine.connect()

        ## gene2go input
        if type(gene2go) == type({}):
            self.gene2go = gene2go
        elif os.path.exists(gene2go):
            tmp = open(gene2go,'r')
            self.gene2go,self.go2gene = cPickle.load(tmp)
            tmp.close()
        else:
            raise Exception("Argument 'gene2go' must be dictionary or pickle file path for (gene2go,go2gene)")

        ## geneset file input
        if not os.path.exists(genesetFile):
            raise Exception("Argument 'genesetFile' must be a valid file path")

        self.read_geneset_file(genesetFile)

    def read_geneset_file(self,genesetFile):
        """
        read in geneset file produced by htsint
        """

        fid = open(genesetFile,'r')
        reader = csv.reader(fid)
        header = reader.next()

        geneset2gene = {}

        if len(header) > 2:
            geneset2transcript = {}
            gene2transcript = {}
            hasTranscripts = True
        else:
            geneset2transcript = None
            gene2transcript = None
            hasTranscripts = False

        for linja in reader:

            ## handle new genesets
            if not geneset2gene.has_key(linja[0]):
                geneset2gene[linja[0]] = set([])
                if hasTranscripts:
                    geneset2transcript[linja[0]] = set([])
            
            ## handle new genes
            if hasTranscripts and not gene2transcript.has_key(linja[1]):
                gene2transcript[linja[1]] = set([])

            ## update the dictionaries
            geneset2gene[linja[0]].update([linja[1]])

            if hasTranscripts:
                geneset2transcript[linja[0]].update(linja[2].split(";"))
            if hasTranscripts:
                gene2transcript[linja[1]].update(linja[2].split(";"))

        ## store as lists
        for key,item in geneset2gene.iteritems():
            geneset2gene[key] = list(item)
        self.geneset2gene = geneset2gene

        if geneset2transcript:
            for key,item in geneset2transcript.iteritems():
                geneset2transcript[key] = list(item)
            self.geneset2transcript = geneset2transcript
        if gene2transcript:
            for key,item in gene2transcript.iteritems():
                gene2transcript[key] = list(item)
            self.gene2transcript = gene2transcript

        fid.close()
        print("-----------------")
        print("%s gene sets loaded from gene set file."%(len(self.geneset2gene.keys())))

    def get_gene_info(self,geneList):
        print("...getting gene info")
        geneInfo = {}
        results = self.conn.execute(Gene.__table__.select(Gene.ncbi_id.in_(geneList)))
        for row in results:
            geneInfo[str(row.ncbi_id)] = {'symbol': str(row.symbol),
                                          'description': str(row.description),
                                          'taxa': self.session.query(Taxon).filter_by(id=row.taxa_id).first().ncbi_id}        
        return geneInfo

    def get_go_terms(self,gene):
        """
        return a human friendly summary of the go terms
        """

        if self.gene2go.has_key(gene):
            terms = [self.session.query(GoTerm).filter(GoTerm.go_id == key).first().name + " (%s)"%key for key in self.gene2go[gene]]
            return ";".join(terms)
        else:
            return 'None'


    def draw_graph(self,geneSetId,layout='spring_layout'):
        """
        create a NetworkX graph 
        """
        
        allLayouts = ['circular_layout',
                      'random_layout',
                      'shell_layout',
                      'spring_layout',
                      'spectral_layout',
                      'fruchterman_reingold_layout']

        if not self.geneset2gene.has_key(geneSetId):
            print("Cannot draw graph for '%s' because ID was not found"%(geneSetId))
            return
        
        print('...fetching information about gene set')
        geneList = self.geneset2gene[geneSetId]
        geneInfo = self.get_gene_info(geneList)
        geneSymbols = []
        for gene in geneList:
            geneSymbols.append(geneInfo[gene]['symbol'])

        ## initialize         
        self.G = nx.Graph()
        for node in geneSymbols:
            self.G.add_node(node)

        ## draw
        print('...drawing and saving')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        pos=nx.spring_layout(self.G)
        nx.draw_networkx_labels(self.G,pos,label_pos=0.3) 
        nx.draw(self.G,ax=ax)
        plt.savefig('foo.png',bbox_inches='tight')

        #nx.draw(Gnode_size=nodeSize1,node_color=colors[0],edge_color=colorList,width=2,edge_cmap=plt.cm.binary,
        #edge_vmin=0.0,edge_vmax=1.0,with_labels=False,alpha=alpha,ax=ax3)



if __name__ == "__main__":
    print "Running..."
