#!/usr/bin/env python
"""
GSA demo script from htsint documentation
"""

import time,os,re
from htsint import GeneOntology,TermDistances,GeneDistances
from htsint.stats import SpectralClustering, SpectralClusterParamSearch
from htsint.blast import BlastMapper

# specify main variables (biological_process, molecular_function, cellular_component)
homeDir = os.path.join(".","demo")
aspect = 'molecular_function' 

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

# Spectral clustering parameter estimation [optional]
#paramSearchPath1 = re.sub("\.csv","-scparams-sv.csv",geneDistancePath)
#paramSearchPath2 = re.sub("\.csv","-scparams-cl.csv",geneDistancePath)
#if not os.path.exists(paramSearchPath1):
#    scps = SpectralClusterParamSearch(geneDistancePath,dtype='distance',aspect=aspect)
#    scps.run(chunks=15)

# Spectral clustering


# Save gene sets



#bm = BlastMapper()
#bmap = bm.load_summary('blast-parsed-summary.csv',best=False)



print("process complete.")
