#!/usr/bin/env python
"""
Demonstrate how to break a pathway into functional modules

"""

import os,sys,re,cPickle
import numpy as np
import xml.sax
from htsint import GeneOntology,TermDistances,GeneDistances
from htsint.stats import SpectralClusterParamSearch,SpectralClusterResults,SpectralCluster

## get a list of the pathways ids
#pathwayList = []
#ORGANISM = "hsa"
#pathways = requests.get('http://rest.kegg.jp/list/pathway/' + ORGANISM)
#for line in pathways.content.split('\n'):
#    pathwayId = line.split('\t')[0].replace('path:', '')
#    pathwayList.append(pathwayId)

## download the txt and xml versions of the pathway
#pathway = "hsa00860"
#ktxt = requests.get('http://rest.kegg.jp/get/' + pathway)
#f = open(os.path.join(".",pathway + ".txt"), 'w')
#f.write(ktxt.content)
#f.close()
#kgml = requests.get('http://rest.kegg.jp/get/' + pathway + '/kgml')
#f = open(os.path.join(".",pathway + '.xml'), 'w')
#f.write(kgml.content)
#f.close()
#print("%s downloaded"%(pathway))

## parse the KEGG pathways
def get_genes(pathway):
    """
    extract genes from a KEGG file
    """
    fid = open(os.path.join(".",pathway+".txt"),'r')
    isGene = False
    gene2symbol = {}
    for linja in fid:
        linja = linja[:-1]
        if re.search("^GENE",linja):
            isGene = True
        elif re.search("^[A-Z]+",linja):
            isGene = False
            
        if isGene:
            geneIds = linja.split(";")[0]
            _geneIds = [re.sub("\s+","",gi) for gi in geneIds.split("  ")]
            geneIds = []
            for gid in _geneIds:
                if len(gid) > 0 and gid != 'GENE':
                    geneIds.append(gid)
            if len(geneIds) != 2:
                print("...%s"%str(geneIds))
                raise Exception("Could not parse gene names")

            ncbiId,symbol = geneIds
            gene2symbol[ncbiId] = symbol
            
    fid.close()
    return gene2symbol

## extract the genes from a KEGG file
pathwayFile = "hsa00860.txt"
pathway = pathwayFile[:-4]
geneList = get_genes(pathway)
print(geneList.items()[:3])

## create a directory for the analysis
gsaDir = os.path.join(".","gsa-path")
if not os.path.exists(gsaDir):
    os.mkdir(gsaDir)

## make imports and specify variables
from htsint import GeneOntology,TermDistances
useIea = True
aspect = "biological_process"
_aspect = 'bp'
taxaList = ['9606']
go = GeneOntology(taxaList,useIea=useIea,aspect=aspect)
termsPath = os.path.join(gsaDir,"go-terms-%s.pickle"%(_aspect))
graphPath = os.path.join(gsaDir,"go-graph-%s.pickle"%(_aspect))

geneIds = geneList.keys()
if not os.path.exists(termsPath):
    go.create_dicts(termsPath,accepted=geneIds)
gene2go,go2gene = go.load_dicts(termsPath)
print("pathway genes with terms: %s/%s"%(len(gene2go.keys()),len(geneIds)))

if not os.path.exists(graphPath):        
    G = go.create_gograph(termsPath=termsPath,graphPath=graphPath)
    print("Term graph for with %s nodes successfully created."%(len(G.nodes())))

# Calculate term distances
termDistancePath = os.path.join(gsaDir,"term-distances-%s.npy"%(_aspect))
if not os.path.exists(termDistancePath):
    td = TermDistances(termsPath,graphPath)
    print("total distances to evaluate: %s"%td.totalDistances)
    td.run_with_multiprocessing(termDistancePath,cpus=7)

# Calculate gene distances
geneDistancePath = os.path.join(gsaDir,"gene-distances-%s.csv"%(_aspect))
if not os.path.exists(geneDistancePath):
    gd = GeneDistances(termsPath,termDistancePath,outFile=geneDistancePath)
    gd.run()

# Spectral Clustering parameter search
silvalFile = re.sub("\.csv","-scparams-sv.csv",geneDistancePath)
clustersFile = re.sub("\.csv","-scparams-cl.csv",geneDistancePath)
if not os.path.exists(silvalFile):
    scps = SpectralClusterParamSearch(geneDistancePath,dtype='distance')
    scps.run(chunks=5,kRange=range(3,11))

## plot the parameter search
psFigureFile = os.path.join(gsaDir,"param-scan-%s.png"%(_aspect))
if not os.path.exists(psFigureFile):
    scr = SpectralClusterResults(silvalFile,clustersFile)
    scr.plot(figName=psFigureFile)

## run spectral clustering
k = 3
sigma = 0.43

labelsPath = os.path.join(gsaDir,"sc-labels-%s.csv"%(_aspect))


if not os.path.exists(labelsPath):
    sc = SpectralCluster(geneDistancePath,dtype='distance')
    sc.run(k,sk=None,sigma=sigma,verbose=True)
    sc.save(labelsPath=labelsPath)


import networkx
from parse_KGML import KGML2Graph
from KeggPathway import KeggPathway
    
p = KeggPathway()
#p.add_node('gene1', data={'type': 'gene', })
#p.get_node('gene1')
#{'type': 'gene'}

graphfile = "%s.xml"%pathway
graph = KGML2Graph(graphfile)[1]

#To get a list of the nodes, use the .nodes() method
print graph.nodes()
print len(graph.nodes())
#print [graph.get_node(n)['label'] for n in graph.nodes()][0:5]
#    ['MGAT1', 'MGAT2', 'C01246', 'TITLE:N-Glycan biosynthesis', 'C03862']
#>>> graph.edges()[0:5]
#    [('56', '57'), ('54', '55'), ('54', '37'), ('54', '58'), ('60', '62')]

#To get detailed informations on a node, use .get_node:
#graph.get_node('10')
#    {'xy': (580, 317), 'type': 'gene', 'label': 'ALG12'}

#All the annotations (such as node type, etc..), are stored in the .label attribute
#graph.label['10']     #doctest: +ELLIPSIS
#{'xy': (580, 317), 'type': 'gene', 'label': 'ALG12'}

#To obtain a subgraph with only the genes of the pathway, it is recommended to use get_genes: 
#genes_graph = graph.get_genes()
#genes_graph.edges()[0:4]
#[('60', '62'), ('63', '3'), ('63', '2'), ('63', '72')]

#for (node1, node2) in genes_graph.edges()[0:4]:
#print genes_graph.get_node(node1)['label'], genes_graph.get_node(node2)
sys.exit()
    
#from xml.dom import minidom

## open xml file
#from xml.dom.minidom import parse
#import xml.dom.minidom

# Open XML document using minidom parser
#DOMTree = xml.dom.minidom.parse("%s.xml"%pathway)
#collection = DOMTree.documentElement

#print dir(collection)
#print collection._get_attributes()
#
#entries = collection.getElementsByTagName("pathway name")
#for entry in entries:
#    print entry
#
#print dir(entry)
