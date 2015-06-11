"""
enrichment functions 
"""

import sys
import numpy as np
from scipy.stats import hypergeom
from htsint.database import fetch_taxa_annotations, db_connect, Taxon, Gene


def enrichment_hypergeo(termList, entityList, species, useIea=True, asGenes=True, aspect='biological_process',verbose=True):

    '''
    termList -- are the terms to be tested
    species  -- an ncbi taxa id
    entityList -- gene or uniprot ids

    What is the probability of finding a given number of terms if we randomly select N out of M objects?

    M -- genes with at least one annotation
    N -- number of draws or size of gene list
    k -- the number of genes annotated by a given term (total type I objects)
    x -- number of times we observe a term in the gene list (draws)

    in R the cdf can be obtained with
    phyper(x,k,M-k,N)
    hypergeom.pmf(x, M, k, N)
    
    Returns a dict where term id is the key and hypergeo pvalue is the value
    '''

    ## connect to db and get annotations for the species
    session, engine = db_connect()
    geneAnnots,uniprotAnnots = fetch_taxa_annotations([species],engine,useIea=useIea,verbose=verbose,aspect=aspect)

    if asGenes == True:
        entity2go = geneAnnots
    else:
        entity2go = uniprotAnnots

    go2entity = {}
    for entity,go in entity2go.iteritems():
        for term in go:
            if not go2entity.has_key(term):
                go2entity[term] = set([])
            go2entity[term].update([entity])
    for go,entity in go2entity.iteritems():
        go2entity[go] = list(entity)

    print('total go terms - %s'%(len(go2entity.keys())))
    print('total entities - %s'%(len(entity2go.keys())))

    ## set variables
    M = len(entity2go.keys())   
    N = len(entityList)         
    results = {}

    for testTerm in termList:
        ## find 
        k = len(go2entity[testTerm])
        x = 0
        for entity in entityList:
            if entity in entity2go and testTerm in entity2go[entity]:
                x += 1

        ## get a p-value
        if 0 in [x,M,N,k]:
            pvalue = np.nan
        else:
            cdf = hypergeom.cdf(x, M, k, N, loc=0)
            if cdf > 0:
                pvalue = 2 * (1-hypergeom.cdf(x, M, k, N))
            else:
                pvalue = 2 * hypergeom.cdf(x, M, k, N)
        results[testTerm] = pvalue

    return results


if __name__ == "__main__":
    print "Running..."

    ## terms to be tested
    termList = ['GO:1902600','GO:0022904','GO:0055114','GO:0042773']

    ## gene list (i.e. differentially expressed cluster)
    geneList = ["13080321","13080320","13080323","13080325","13080324","13080326"]

    species = '5476'
    results = enrichment_hypergeo(termList, geneList, species, useIea=True)
