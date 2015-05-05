#!/usr/bin/env python
"""
A class to represent and draw gene sets
"""

__author__ = "Adam Richards"

import os,cPickle,csv,re
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from htsint.database import db_connect,Gene,Taxon,GoTerm
from htsint.blast import BlastMapper

class GeneSet(object):
    """
    gene set class
    """

    def __init__(self,verbose=True):
        """
        Constructor
        """

        ## setup db
        self.session, self.engine = db_connect(verbose=verbose)
        self.conn = self.engine.connect()

        ## global variables
        self.verbose = verbose
        self.dpi = 400
        self.labelOffset = 0.05
        self.fontSize = 8 
        self.legendFontSize = 8
        self.fontName = 'serif'
        self.termAlpha = 0.95
        self.geneAlpha = 0.4
        self.nodeSizeGene = 300
        self.nodeSizeTerm = 300
        self.colors = ['#000000','#FFCC33','#3333DD']
        self.cmap = plt.cm.Blues
        self.lineEnd = 53
        self.linesMax = 32


    def load_geneset(self,genesetFile,gene2go,distMat):
        """
        gene2go - dictionary or file path
        """

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

        ## distance matrix
        if type(distMat) != type(np.array([])):
            raise Exception("Argument 'distMat' must be a NumPy array")
        self.distMat = distMat

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
        if self.verbose:
            print("-----------------")
            print("%s gene sets loaded from gene set file."%(len(self.geneset2gene.keys())))

    def get_gene_info(self,geneList):
        if self.verbose:
            print("...getting gene info")
        geneInfo = {}
        results = self.conn.execute(Gene.__table__.select(Gene.ncbi_id.in_(geneList)))
        for row in results:
            taxaQuery = self.session.query(Taxon).filter_by(id=row.taxa_id).first()
            geneInfo[str(row.ncbi_id)] = {'symbol': str(row.symbol),
                                          'description': str(row.description),
                                          'taxa': taxaQuery.ncbi_id,
                                          'species': taxaQuery.name
                                      }
        return geneInfo

    def get_term_distances(self,termList,percentile=50):
        """
        return a dictionary of term distances
        ignores any distances greater than percentile threshold
        """

        def scale(val, src, dst):
            """ 
            Scale the given value from the scale of src to the scale of dst. 
            """
            return ((val - src[0]) / (src[1]-src[0])) * (dst[1]-dst[0]) + dst[0]

        mat = self.distMat[:,2].astype(float)

        ## get the percentile threshold
        threshold = np.percentile(mat,percentile)
        if self.verbose:
            print "percentile threshold: %s (%s)"%(threshold,percentile)

        termDist = {}
        for i in range(self.distMat.shape[0]):
            linja = self.distMat[i,:]
            if linja[0] not in termList or linja[1] not in termList:
                continue
            if float(linja[2]) > threshold:
                continue
            if not termDist.has_key(linja[0]):
                termDist[linja[0]] = {}
            if not termDist[linja[0]].has_key(linja[1]):
                termDist[linja[0]][linja[1]] = 1.0 - scale(float(linja[2]),(mat.min(),mat.max()),(0,1))

        return termDist

    def get_go_term(self,term):
        """
        return a human friendly summary of the go terms
        """

        description =self.session.query(GoTerm).filter(GoTerm.go_id == term).first().name
        if description:
            return description
        else:
            return 'None'

    def get_term_counts(self,geneList):
        
        ## count the terms
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
        termIds = np.array(["%s"%str(i+1) for i in range(termList.size)])

        toReturn = []
        for r,rank in enumerate(rankedInds):
            toReturn.append({'rank':r+1,'term':termList[rank],'desc':self.get_go_term(termList[rank]),'counts':edgeCounts[rank]})

        return toReturn

    def draw_network(self,geneSetId,layout='spring',name='None',percentile=50,ax=None):
        """
        create a NetworkX graph 

        allLayouts = ['circular_layout',
                      'random_layout',
                      'shell_layout',
                      'spring_layout',
                      'spectral_layout',
                      'fruchterman_reingold_layout']
        """
                
        if not name:
            pass
        if name == 'None':
            name = "%s_network.png"%(geneSetId)

        if not self.geneset2gene.has_key(geneSetId):
            print("Cannot draw graph for '%s' because ID was not found"%(geneSetId))
            return
        
        if self.verbose:
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
        termIds = np.array(["%s"%str(i+1) for i in range(termList.size)])

        term2id = {}
        for r,rank in enumerate(rankedInds):
            term2id[termList[rank]] = termIds[r]

        termIds = termIds.tolist()    
        geneInfo = self.get_gene_info(geneList)
        taxonNames = list(set([geneInfo[gene]['taxa'] for gene in geneList]))
        taxonNames.sort()
        taxonIds = {}
        for t,taxon in enumerate(taxonNames):
            if self.verbose:
                print("%s,%s"%(t,taxon))
            taxonIds[taxon] = str(t)

        geneSymbols = {}
        gene2symbol = {}
        geneLabels = {}
        for gene in geneList:
            nodeId = geneInfo[gene]['symbol'] + "-" + taxonIds[geneInfo[gene]['taxa']]
            taxonId = geneInfo[gene]['taxa']
            gene2symbol[gene] = nodeId
            if not geneLabels.has_key(taxonId):
                geneSymbols[taxonId] = []
                geneLabels[taxonId] = {}
            geneLabels[taxonId][nodeId] = geneInfo[gene]['symbol']
            geneSymbols[taxonId].append(nodeId)

        ## annotation edges
        edgeList1 = [(gene2symbol[edge[0]],term2id[edge[1]]) for edge in geneEdges]
        
        ## term term edges
        termDistances = self.get_term_distances(termList,percentile=percentile)

        edgeList2 = []
        for term1,items in termDistances.iteritems():
            for term2,dist in items.iteritems():
                edgeList2.append((term2id[term1],term2id[term2],dist))

        ## initialize
        self.G = nx.Graph()
        for taxon in taxonNames:
            for node in geneSymbols[taxon]:
                self.G.add_node(node)
        for term in termIds:
            self.G.add_node(term)
        for edge in edgeList1:
            self.G.add_edge(edge[0],edge[1],weight=3.0)
        edge2colors = []
        for edge in edgeList2:
            self.G.add_edge(edge[0],edge[1],weight=edge[2])
            edge2colors.append(edge[2])
        if self.verbose:
            print("...there are %s annotation edges"%(len(edgeList1)))
            print("...there are %s term-term edges"%(len(edgeList2)))
            print('...drawing and saving')

        ## use provided axis or create a new one
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_yticks([])
            ax.set_xticks([])

        ## layout
        if layout == 'spring':
            pos = nx.spring_layout(self.G,scale=1)
        elif layout == 'spectral':
            pos = nx.spectral_layout(self.G)
        else:
            print("layout not supported %s... using spring"%(layout))
            pos = nx.spring_layout(self.G)

        ## draw gene nodes for each species
        colorItr = 1

        for taxon in taxonNames:
            symbols = geneSymbols[taxon]
            nodeSizes = [self.nodeSizeGene + (15 * float(len(s))) for s in symbols]
            _G = self.G.subgraph(symbols)
            _pos = {}
            for p in pos:
                if p in symbols:
                    _pos[p] = pos[p]

            nx.draw_networkx_nodes(_G,_pos,node_size=nodeSizes,nodelist=symbols,node_shape=r'o',
                                   node_color=self.colors[colorItr],alpha=self.geneAlpha,ax=ax)
            colorItr += 1

            ## offset labels
            #for p in _pos:
            #    _pos[p][1] += self.labelOffset
            nx.draw_networkx_labels(_G,_pos,font_color='black',font_family=self.fontName,labels=geneLabels[taxon],ax=ax)

        ## draw term nodes
        nx.draw_networkx_nodes(self.G,pos,node_size=self.nodeSizeTerm,nodelist=termIds,node_shape='s',
                               node_color=self.colors[0],alpha=self.termAlpha,ax=ax)
        ## term labels
        G1 = self.G.subgraph(termIds)
        pos1 = {}
        for p in pos:
            if p in termIds:
                pos1[p] = pos[p]
        nx.draw_networkx_labels(G1,pos1,font_color='white',font_family=self.fontName,font_size=self.fontSize,ax=ax)

        ## draw edges
        nx.draw_networkx_edges(self.G,pos,edgelist=edgeList1,width=1.0,edge_color='k',style='dashed',ax=ax,alpha=0.3)
        nx.draw_networkx_edges(self.G,pos,edgelist=edgeList2,edge_color=edge2colors,width=2.0,style='solid',edge_cmap=self.cmap,ax=ax)

        if name:
            plt.savefig(name,bbox_inches='tight',dpi=self.dpi)

        return term2id,termList,taxonNames

    def draw_figure(self,geneSetId,layout='spring',name=None,percentile=50):

        ## specify axes l,b,w,h 
        fig = plt.figure(figsize=(10,6))
        ax1  = fig.add_axes([0.3, 0.0, 0.7, 1.0])
        ax2  = fig.add_axes([0.0, 0.07, 0.3, 0.93],axisbg='#EEEEEE')
        ax3  = fig.add_axes([0.0, 0.0, 0.3, 0.07],axisbg='#EEEEEE')      # empty
        ax4  = fig.add_axes([0.015, 0.03, 0.27, 0.03])

        for ax in [ax1,ax2,ax3,ax4]:
            ax.set_yticks([])
            ax.set_xticks([])
        
        ax2.set_ylim((0,1))
        ax2.set_xlim((0,1))
        #ax1.set_frame_on(False)

        def add_line(toPrint,lineCount,current):
            toPrint = toPrint[:self.lineEnd]
            if lineCount == self.linesMax:
                ax2.text(0.01,current," ...",color='k',fontsize=self.legendFontSize,fontname=self.fontName,ha="left", va="center")
                current = current - increment
                lineCount += 1
                return current, lineCount
            elif lineCount > self.linesMax:
                return current, lineCount

            ax2.text(0.01,current,toPrint,color='k',fontsize=self.legendFontSize,fontname=self.fontName,ha="left", va="center")
            current = current - increment
            lineCount += 1
        
            return current, lineCount

        ## draw the network
        term2id,termList,taxonNames = self.draw_network(geneSetId,layout=layout,name=None,percentile=percentile,ax=ax1)

        ## axis 2 (legend)
        lineCount = 0
        current = 0.98
        increment = 0.03
        G1 = nx.Graph()
        colorIter = 1
        pos1 = {}
        for node in taxonNames:
            G1.add_node(node)
            pos1[node] = (0.05,current)
            nx.draw_networkx_nodes(G1,pos1,node_size=100,nodelist=[node],
                                   node_color=self.colors[colorIter],alpha=self.geneAlpha,ax=ax2)
            query = self.session.query(Taxon).filter_by(ncbi_id=node).first()
            toAdd = "        %s (%s)"%(query.name,node)[:self.lineEnd]
            current,lineCount = add_line(toAdd,lineCount,current)
            colorIter += 1
        
        current,lineCount = add_line("-"*self.lineEnd,lineCount,current)
        termIds = [int(i) for i in list(set(term2id.values()))]
        termIds.sort()
        termIds = [str(i) for i in termIds]
        id2term = {}
        for key,item in term2id.iteritems():
            id2term[item] = key
        
        for termId in termIds:
            term = id2term[termId]
            desc = self.get_go_term(term)
            toPrint = "%s - %s"%(termId.zfill(3),desc)
            #print toPrint, term
            
            ## check to see if we need to use multiple lines
            if len(toPrint) > self.lineEnd:
                wordBreaks = np.array([m.start(0) for m in re.finditer("\s+",toPrint)])
                lineBreak = wordBreaks[np.where(wordBreaks < self.lineEnd)[0]][-1]
                remaining = toPrint[lineBreak:]
                toPrint = toPrint[:lineBreak]
            else:
                remaining = None
            toPrint = toPrint[:self.lineEnd]
            current,lineCount = add_line(toPrint,lineCount,current)

            if remaining:
                current,lineCount = add_line("        "+remaining,lineCount,current)

        ## axis 4 (colorbar)
        norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
        cbar = mpl.colorbar.ColorbarBase(ax4,cmap=self.cmap,
                                         ticks=[0.0,0.25,0.5,0.75,1.0],
                                         norm=norm,
                                         orientation='horizontal')
        cbar.ax.tick_params(labelsize=self.fontSize) 

        ## save
        plt.savefig(name,bbox_inches='tight',dpi=self.dpi)

if __name__ == "__main__":
    print "Running..."
