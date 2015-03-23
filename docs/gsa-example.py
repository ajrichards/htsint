#!/usr/bin/env python
"""
GSA demo script from htsint documentation
"""

import time
from htsint.blast import BlastMapper
from htsint import GeneOntology,TermDistances,GeneDistances

bm = BlastMapper()
bmap = bm.load_summary('blast-parsed-summary.csv',best=False)

# Create a term graph (biological_process, molecular_function, cellular_component)
go = GeneOntology(["8364","8355"],useIea=False,aspect='molecular_function')
termsPath = "go-terms.pickle"
graphPath = "go-graph.pickle"
go.create_dicts(termsPath)
gene2go,go2gene = go.load_dicts(termsPath)
G = go.create_gograph(termsPath=termsPath,graphPath=graphPath)
print("%s genes have at least one annotation"%(len(gene2go.keys())))
print("Term graph for with %s nodes successfully created."%(len(G.nodes())))

# Calculate term distances
td = TermDistances(termsPath,graphPath)
print("total distances to evaluate: %s"%td.totalDistances)
timeStart = time.time()
termDistancePath = "term-distances.npy"
td.run_with_multiprocessing(termDistancePath,cpus=16)
print time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))

# Calculate gene distances
#geneDistancePath = "gene-distances.csv"
#gd = GeneDistances(termsPath,graphPath,termDistancePath,outFile=geneDistancePath)
#gd.run()

# Spectral Clustering



# Save gene sets

