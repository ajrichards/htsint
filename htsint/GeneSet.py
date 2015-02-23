#!/usr/bin/env python
"""
A class to represent and draw gene sets
"""

__author__ = "Adam Richards"

import os,cPickle,csv
import numpy as np
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

    def draw_graph(self,geneSetId,layout='spring_layout',name=None):
        """
        create a NetworkX graph 

        allLayouts = ['circular_layout',
                      'random_layout',
                      'shell_layout',
                      'spring_layout',
                      'spectral_layout',
                      'fruchterman_reingold_layout']
        """
        
        nodeSizeGene = 500
        nodeSizeTranscript = 100
        nodeSizeTerm = 200
        colors = ['#FFCC33','#3333DD','#000000']

        if not name:
            name = "%s_network.png"%(geneSetId)

        if not self.geneset2gene.has_key(geneSetId):
            print("Cannot draw graph for '%s' because ID was not found"%(geneSetId))
            return
        
        print('...fetching information about gene set')
        geneList = self.geneset2gene[geneSetId]

        ## get terms
        geneEdges = []
        termCounts = {}
        for gene in geneList:
            if not self.gene2go.has_key(gene):
                continue
            for term in self.gene2go[gene]:
                geneEdges.append((gene,term))
                if not termCounts.has_key(term):
                    termCounts[term] = 0
                termCounts[term] += 1

        ## rank the terms by the number of connections
        termList = np.array(termCounts.keys())
        edgeCounts = np.array(termCounts.values())
        rankedInds = np.argsort(edgeCounts)[::-1]
        termIds = np.array(["GO:%s"%str(i+1) for i in range(termList.size)])
        print termIds

        term2id = {}
        for r,rank in enumerate(rankedInds):
            print termIds[r],termList[rank],edgeCounts[rank]
            term2id[termList[rank]] = termIds[r]
            
        geneInfo = self.get_gene_info(geneList)
        geneSymbols = []
        gene2symbol = {}
        for gene in geneList:
            gene2symbol[gene] = geneInfo[gene]['symbol']
            geneSymbols.append(geneInfo[gene]['symbol'])
   
        edgeList = [(gene2symbol[edge[0]],term2id[edge[1]]) for edge in geneEdges]

        ## initialize
        self.G = nx.Graph()
        for node in geneSymbols:
            self.G.add_node(node)
        for term in termIds:
            self.G.add_node(term)
        for edge in edgeList:
            self.G.add_edge(edge[0],edge[1],weight=1)
    
        ## draw
        print('...drawing and saving')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        if layout == 'spring':
            pos1 = nx.spring_layout(self.G)
        else:
            print("layout not supported %s... using spring"%(layout))
            pos1 = nx.spring_layout(self.G)
        
        print termList

        nx.draw_networkx_nodes(self.G,pos1,node_size=nodeSizeGene,nodelist=geneSymbols,node_shape='',
                               node_color=colors[0],ax=ax)
        nx.draw_networkx_nodes(self.G,pos1,node_size=nodeSizeTerm,nodelist=termIds.tolist(),node_shape='s',
                               node_color=colors[1],font_color='white',ax=ax)
        nx.draw_networkx_labels(self.G,pos1,nodelist=geneSymbols,label_pos=5.0,ax=ax,font_color='black')
        nx.draw_networkx_labels(self.G,pos1,nodelist=termIds.tolist(),label_pos=5.0,ax=ax,font_color='white')
        nx.draw_networkx_edges(self.G,pos1,edgelist=edgeList,width=1,edge_color='k',style='solid',ax=ax)

        plt.savefig(name,bbox_inches='tight',dpi=400)


if __name__ == "__main__":
    print "Running..."
