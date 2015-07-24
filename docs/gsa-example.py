#!/usr/bin/env python
"""
GSA demo script from htsint documentation
"""

import time,os,re
from htsint import GeneOntology,TermDistances,GeneDistances,GeneSetCollection
from htsint.stats import SpectralCluster, SpectralClusterParamSearch, SpectralClusterResults
from htsint.blast import BlastMapper

# specify main variables (biological_process, molecular_function, cellular_component)
homeDir = os.path.join(".","demo")
if not os.path.isdir(homeDir):
    os.mkdir(homeDir)

_aspect = 'bp'
aspect = 'biological_process' 

# Create a term graph
go = GeneOntology(["8364","8355"],useIea=False,aspect=aspect)
termsPath = os.path.join(homeDir,"go-terms.pickle")
graphPath = os.path.join(homeDir,"go-graph.pickle")
if not os.path.exists(termsPath):
    go.create_dicts(termsPath)
gene2go,go2gene = go.load_dicts(termsPath)
G = go.create_gograph(termsPath=termsPath,graphPath=graphPath)
print("%s genes have at least one annotation"%(len(gene2go.keys())))
print("Term graph for with %s nodes successfully created."%(len(G.nodes())))

# Calculate term distances (most time consuming step)
termDistancePath = os.path.join(homeDir,"term-distances.npy")
if not os.path.exists(termDistancePath):
    td = TermDistances(termsPath,graphPath)
    print("total distances to evaluate: %s"%td.totalDistances)
    timeStart = time.time()
    td.run_with_multiprocessing(termDistancePath,cpus=30)
    print time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))

# Calculate gene distances
geneDistancePath = os.path.join(homeDir,"gene-distances.csv")
if not os.path.exists(geneDistancePath):
    gd = GeneDistances(termsPath,termDistancePath,outFile=geneDistancePath)
    gd.run()

# Spectral Clustering parameter search 
silvalFile = re.sub("\.csv","-scparams-sv.csv",geneDistancePath)
clustersFile = re.sub("\.csv","-scparams-cl.csv",geneDistancePath)
if not os.path.exists(silvalFile):
    scps = SpectralClusterParamSearch(geneDistancePath,dtype='distance')
    scps.run(chunks=15)

## plot the parameter search 
psFigureFile = os.path.join(homeDir,"param-scan-%s.png"%(_aspect))
if not os.path.exists(psFigureFile):
    scr = SpectralClusterResults(silvalFile,clustersFile)
    scr.plot(figName=psFigureFile)

## run spectral clustering
k = 20
sigma = 0.08

labelsPath = os.path.join(homeDir,"sc-labels-%s.csv"%(_aspect))
if not os.path.exists(labelsPath):
    sc = SpectralCluster(geneDistancePath,dtype='distance')
    sc.run(k,sk=None,sigma=sigma,verbose=True)
    sc.save(labelsPath=labelsPath)

## Save gene sets
bm = BlastMapper()
bmap = bm.load_summary('blast-parsed-summary.csv',best=False,taxaList=['8355','8364'])

transcriptMin,transcriptMax = 9,1000  
gsFile = os.path.join(homeDir,"%s.gmt"%(_aspect))                                                                                                       
if not os.path.exists(gsFile):
    gsc = GeneSetCollection(labelsPath,gene2go)
    gsc.write(blastMap=bmap,transcriptMin=transcriptMin,transcriptMax=transcriptMax,outFile=gsFile)

print("process complete.")
