#!/usr/bin/env python
"""
Demonstrate how to break a pathway into functional modules

"""

import os,sys,re,cPickle
import numpy as np
from htsint import GeneOntology,TermDistances

## import requests
# down load all ecoli pathways in xml
#ORGANISM = "hsa"
#pathwayList = []
#pathways = requests.get('http://rest.kegg.jp/list/pathway/' + ORGANISM)
#for line in pathways.content.split('\n'):
#    pathwayId = line.split('\t')[0].replace('path:', '')
#    pathwayList.append(pathwayId)
#    ktxt = requests.get('http://rest.kegg.jp/get/' + pathwayId)
#    f = open(os.path.join("pathways",pathwayId + ".txt"), 'w')
#    f.write(ktxt.content)
#    #kgml = requests.get('http://rest.kegg.jp/get/' + pathwayId + '/kgml')
#    #f = open(os.path.join("pathways",pathwayId + '.xml'), 'w')
#    #f.write(kgml.content)
#    f.close()

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
print geneList

## create a directory for the analysis
gsaDir = os.path.join(".","gsa")
if not os.path.exists(gsaDir):
    os.mkdir(gsaDir)

sys.exit()

## CREATE FUNCTIONAL MODULES FOR PATHWAY (LOOP OVER SPECIES OF INTEREST)
# Create a term graph (biological_process, molecular_function, cellular_component)

useIea = True
aspect = "biological_process"
_aspect = 'bp'
taxaList = ['9606']
go = GeneOntology(taxaList,useIea=useIea,aspect=aspect)
termsPath = os.path.join(gsaDir,"go-terms-%s-%s.pickle"%(_aspect))
graphPath = os.path.join(gsaDir,"go-graph-%s-%s.pickle"%(_aspect))

if not os.path.exists(termsPath):
    go.create_dicts(termsPath)
    _gene2go,_go2gene = go.load_dicts(termsPath)

    actualGeneList = []
    actualGeneClasses = []
    gene2go = {}
    go2gene = {}
    for g,gene in enumerate(geneList):
        if _gene2go.has_key(gene):
            actualGeneList.append(gene)
            actualGeneClasses.append(geneClasses[g])
            terms =  _gene2go[gene]
            gene2go[gene] = terms
            
            for term in terms:
                if go2gene.has_key(term) == False:
                    go2gene[term] = set([])
                go2gene[term].update([gene])

    for term,genes in go2gene.iteritems():
        go2gene[term] = list(genes)

    ## overwrite dictionaries
    tmp = open(termsPath,'w')
    cPickle.dump([gene2go,go2gene],tmp)
    tmp.close()
        
    print "original genes: %s/%s"%(len(geneList), len(_gene2go.keys()))
    print "actual list: %s/%s"%(len(actualGeneList), len(gene2go.keys()))
        
    G = go.create_gograph(termsPath=termsPath,graphPath=graphPath)
    print("%s genes have at least one annotation"%(len(gene2go.keys())))
    print("Term graph for with %s nodes successfully created."%(len(G.nodes())))

# Calculate term distances
termDistancePath = os.path.join(gsaDir,"term-distances-%s-%s.npy"%(expId,aspect))
if not os.path.exists(termDistancePath):
    td = TermDistances(termsPath,graphPath)
    print("total distances to evaluate: %s"%td.totalDistances)
    td.run_with_multiprocessing(termDistancePath,cpus=7)
    
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
                                    
