#!/usr/bin/env python

import os,sys,csv,re
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene


def get_blast_map(resultsFilePath,evalue=0.00001,taxaList=None,asGenes=False,append=False):
    """
    load assembly blast results into dictionary

    if taxaList is provided then only genes from given taxa will be included in map

    if asGene == True the results are provided with keys to genes not isoforms


    """

    if not os.path.exists(resultsFilePath):
        raise Exception("cannot find results file path %s"%resultsFilePath)

    if taxaList != None:
        ## prepare database connections
        session,engine = db_connect()
        conn = engine.connect()
        s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(taxaList))
        _taxaQueries = conn.execute(s)
        taxaQueries = _taxaQueries.fetchall()
        totalQueries = set([])
        filteredQueries = set([])
        filteredHits = set([])
        selectedTaxa = [str(tquery['id']) for tquery in taxaQueries]
        taxa2name = dict([(str(tquery['id']),str(tquery['ncbi_id'])) for tquery in taxaQueries])
        
        ## create a gene2taxa dictionary
        #gene2taxa,gene2desc = {},{}
        #for tquery in taxaQueries:
        #    s = select([Gene.taxa_id,Gene.ncbi_id,Gene.description],Gene.taxa_id==tquery['id'])
        #    _geneQueries = conn.execute(s)
        #    geneQueries = _geneQueries.fetchall()
        #    gene2taxa.update(dict([(str(r['ncbi_id']),str(r['taxa_id'])) for r in geneQueries]))
        #    gene2desc.update(dict([(str(r['ncbi_id']),str(r['description'])) for r in geneQueries]))

    results = {}
    fid = open(resultsFilePath,'rU')
    reader = csv.reader(fid)
    header = reader.next()
    print header

    ## loop through file and save best
    uniqueQueries = set([])
    totalQueries = 0
    unfilteredQueries = 0

    for linja in reader:
        
        if len(linja) == 4:
            queryId = linja[0]
            hitId = linja[1]
            hitNcbiId = linja[2]
            _evalue = float(linja[3])
        else:
            queryId = linja[0]
            queryNcbi = linja[1]
            hitId = linja[2]
            hitNcbiId = linja[3]
            _evalue = linja[4]

        if asGenes == True:
            queryId = re.sub("_i\d+","",queryId)

        # filtering
        totalQueries += 1
        if '-' in linja:
            continue
        if _evalue > evalue:
            continue

        if taxaList and gene2taxa.has_key(str(hitNcbiId)) == False:
            continue

        unfilteredQueries += 1
        uniqueQueries.update([queryId])
        
        ## use the best evalue
        if not results.has_key(queryId):
            if append:
                results[queryId] = [(hitNcbiId,_evalue)]
            else:
                results[queryId] = (hitNcbiId,_evalue)
        if _evalue < results[queryId][1]:
            if append:
                results[queryId].append((hitNcbiId,_evalue))
            else:
                results[queryId] = (hitNcbiId,_evalue)
                
    uniqueQueries = list(uniqueQueries)
    
    print("total queries: %s"%totalQueries)
    print("unfiltered queries: %s"%unfilteredQueries)
    print("unique: %s"%len(uniqueQueries))
        
    return results


def create_blast_map(refTaxon,taxaList,resultsFilePath,evalue=0.00001,verbose=False):
    """
    read a summarized reference blast results file and create a map
    results are gene centric

    example results file looks like this
    query(refseq),query(geneId),hit(uniprotEntry),hit(geneEntry),e-value
    NP_001016845.1,549599,AQP3_HUMAN,360,1.24637e-170
    NP_001016845.1,549599,AQP9_HUMAN,366,7.57313e-92
    NP_001016845.1,549599,AQP10_HUMAN,89872,2.87154e-85
    NP_001016845.1,549599,AQP7_HUMAN,364,8.01308e-84
    NP_001016845.1,549599,AQP5_HUMAN,362,4.267e-15

    """

    ## error check
    if not os.path.exists(resultsFilePath):
        raise Exception("cannot find results file path %s"%resultsFilePath)

    if refTaxon not in taxaList:
        raise Exception("refTaxon must be in taxaList")

    ## prepare database connections
    session,engine = db_connect()
    conn = engine.connect()

    ## read through the file to map the genes to taxa ids    
    s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(taxaList))
    _taxaQueries = conn.execute(s)
    taxaQueries = _taxaQueries.fetchall()
    totalQueries = set([])
    filteredQueries = set([])
    filteredHits = set([])
    selectedTaxa = [str(tquery['id']) for tquery in taxaQueries]
    taxa2name = dict([(str(tquery['id']),str(tquery['ncbi_id'])) for tquery in taxaQueries])

    ## create a gene2taxa dictionary
    gene2taxa = {}
    for tquery in taxaQueries:
        s = select([Gene.taxa_id,Gene.ncbi_id],Gene.taxa_id==tquery['id'])
        _geneQueries = conn.execute(s)
        gene2taxa.update(dict([(str(r['ncbi_id']),str(r['taxa_id'])) for r in _geneQueries.fetchall()]))

    ## creats a dictionary results['geneId]['taxaId'] = bestHitMappedTaxa
    results = {}
    fid = open(resultsFilePath,'rU')
    reader = csv.reader(fid)
    header = reader.next()
        
    ## loop through file and save best
    for linja in reader:
        _evalue = float(linja[4])
        totalQueries.update([linja[1]])

        ## filter by species and evalue
        if not gene2taxa.has_key(linja[3]) or _evalue > evalue:
            continue
        if linja[1] == '-' or linja[3] == '-':
            continue

        ## filter self matches
        if taxa2name[gene2taxa[linja[3]]] == refTaxon:
            continue
        if taxa2name[gene2taxa[linja[1]]] != refTaxon:
            raise Exception("Invalid query or invalid refTaxon %s != %s"%(taxa2name[gene2taxa[linja[1]]],refTaxon))
    
        taxId = taxa2name[gene2taxa[linja[3]]]
        filteredQueries.update([linja[1]])
        filteredHits.update([linja[3]])
        
        if not results.has_key(taxId):
            results[taxId] = {}
        
        ## use the best evalue
        if not results[taxId].has_key(linja[1]):
            results[taxId][linja[1]] = (linja[3],_evalue)
        if _evalue < results[taxId][linja[1]][1]:
            results[taxId][linja[1]] = (linja[3],_evalue)

    fid.close()
    
    ## returns a simplified form of the results as two mappers
    mapper1,mapper2 = {},{}
    test1 = set([])
    for taxId in results.iterkeys(): 
        for queryGene,hit in results[taxId].iteritems():
            if not mapper1.has_key(hit[0]):
                mapper1[hit[0]] = []
            mapper1[hit[0]].append(queryGene)
            
            if not mapper2.has_key(queryGene):
                mapper2[queryGene] = []
            mapper2[queryGene].append(hit[0])

    #debug = set([])
    #for key, item in mapper.iteritems():
    #    debug.update(item)
    
    #for tquery in taxaQueries:
    #    taxId = str(tquery['ncbi_id'])
    #    for key,item in results[taxId].iteritems():
    #        debug.update([key])
    #        mapper[item[0]] = key

    #debug = list(debug)
    #print 'debug', len(debug),missing
    #print results.keys(),len(results['8364'].keys()), len(results['8355'].keys()), len(list(set(results['8364'].keys() + results['8355'].keys())))
    print('BLAST: total queries: %s'%(len(list(totalQueries))))
    print('BLAST: filtered queries (evalue=%s)(taxa=%s): %s'%(evalue,str(taxaList),len(list(filteredQueries))))
    print('BLAST: filtered hits: %s'%(len(list(filteredHits))))
    return mapper1,mapper2
