#!/usr/bin/env python
"""
library of functions for use with the GeneOntology class
"""

import os,sys,re,time
from sqlalchemy.sql import select
from htsint import Configure
from DatabaseTables import Taxon,Uniprot,Gene,GoTerm,GoAnnotation

def remove_empty(lst):
    if None in lst:
        lst.remove(None)

def get_ontology_file():
    """
    check for presence of ontology file
    raise exception when not found
    return the file path
    """

    config = Configure()
    dataDir = config.log['data']
    ontologyFile = os.path.join(dataDir,'go.obo')
    if os.path.exists(ontologyFile) == False:
        raise Exception("Could not find 'go.obo' -- did you run DatabaseFetch?")

    return ontologyFile

def read_ontology_file():
    """
    read the ontology file to find term-term edges
    store all the relationships in a dictionary
    """

    ontologyFile = get_ontology_file()
    fid = open(ontologyFile,'r')
    termCount = 0
    goId = None
    goDict = {"cellular_component":{},
              "molecular_function":{},
              "biological_process":{}}

    def add_term(goNamespace,source,sink):
        if len(sink) != 1 or source == sink[0]:
            return
        if not re.search("GO\:",source) or not re.search("GO\:",sink[0]):
            raise Exception("Invalid go id in ontology file: %s, %s"%(source,sink[0]))

        if goDict[goNamespace].has_key(source) == False:
            goDict[goNamespace][source] = set([])
        goDict[goNamespace][source].update(sink)

    for linja in fid.readlines():
        linja = linja[:-1]

        ## find go id
        if re.search("^id\:",linja):
            goId = re.sub("^id\:|\s+","",linja)
            goName,goNamespace = None,None
            is_a,part_of = None,None
            isObsolete = False
            termCount += 1
            continue 

        ## find a few go id attributes
        if re.search("^name\:",linja):
            goName = re.sub("^name\:\s+","",linja)
        if re.search("^namespace\:",linja):
            goNamespace = re.sub("^namespace\:\s+","",linja)
        if re.search("^def\:",linja):
            goDef = re.sub("^def\:\s+","",linja)
            if re.search("OBSOLETE\.",goDef):
                isObsolete = True

        ## ignore obolete terms
        if goId == None or isObsolete == True:
            continue

        ## add term relationships (is_a)
        if re.search("^is\_a\:",linja):
            is_a = re.sub("^is\_a\:\s+","",linja)
            is_a_sink = re.findall("GO\:\d+",is_a)
            add_term(goNamespace,goId,is_a_sink)

        ## other term relationships
        #if re.search("^relationship\:",linja):
        #    part_of_sink = re.findall("GO\:\d+",linja)
        #    add_term(goNamespace,goId,part_of_sink)
            
    return goDict

def get_annotation_file():
    """
    check for presence of the annotation file
    raise exception when not found
    return the file path
    """

    fileName = 'gene_association.goa_uniprot.db'
    config = Configure()
    dataDir = config.log['data']
    annotationFile = os.path.join(dataDir,fileName)
    if os.path.exists(annotationFile) == False:
        raise Exception("Could not find '%s' -- did you run FetchDbData.py?"%(fileName))

    return annotationFile

def get_total_annotations():
    """
    get the number of annotations in the uniprot file
    this does not include the gene2go file
    """

    config = Configure()
    taxaList = config.log['taxa']

    annotationFile = get_annotation_file()
    annotationFid = open(annotationFile,'rU')
    annotsCount = 0
    annotatedIds = {}
    totalAnnotations = 0

    for record in annotationFid:
        record = record[:-1].split("\t")
        if record[0][0] == "!":
            continue
        if record[0] != 'UniProtKB':
            continue

        taxon = re.sub("taxon:","",record[12])
        if taxon == "" or re.search("\|",taxon):
            continue
        if taxon not in taxaList:
            continue

        totalAnnotations += 1

    gene2goFile = get_gene2go_file()
    gene2goFid = open(gene2goFile,'rU')
    header = gene2goFid.next()

    for record in gene2goFid:
        totalAnnotations += 1

    return totalAnnotations

def get_gene2go_file():
    """
    check for presence of the annotation file
    raise exception when not found
    return the file path
    """

    config = Configure()
    dataDir = config.log['data']
    annotationFile = os.path.join(dataDir,'gene2go.db')
    if os.path.exists(annotationFile) == False:
        raise Exception("Could not find 'gene2go' -- did you run FetchDbData.py?")

    return annotationFile

def get_evidence_codes(useIea=False):
    """
    returns the list of evidence codes to be used by default
    """

    expEvidCodes = ["EXP","IDA","IPI","IMP","IGI","IEP"]
    compEvidCodes = ["ISS","ISO","ISA","ISM","IGC","RCA"]
    statEvidCodes = ["TAS","NAS","IC"]
    nonCuratedEvidCodes = ["IEA"]

    if useIea == False:
        return expEvidCodes + statEvidCodes + compEvidCodes
    else:
        return expEvidCodes + statEvidCodes + compEvidCodes + nonCuratedEvidCodes

def get_annotated_genes(session,taxonId,useIea=False,aspect='biological_process',filterTaxa=True):
    """
    fetch annotations by taxa
    """

    if aspect not in ['biological_process','cellular_component','molecular_function']:
        raise Exception("Invalid aspect specified")

    taxaQuery = session.query(Taxon).filter_by(ncbi_id=taxonId).first()
    acceptedCodes = get_evidence_codes(useIea=useIea)
    annotations = session.query(GoAnnotation).join(GoTerm).\
                  filter(GoAnnotation.taxa_id==taxaQuery.id).\
                  filter(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                  filter(GoTerm.aspect==aspect).all()

    annotatedGenes = list(set([a.gene_id for a in annotations]))
    annotatedProts = list(set([a.uniprot_id for a in annotations]))
    remove_empty(annotatedGenes)
    remove_empty(annotatedProts)
    
    if filterTaxa == False:
        return [session.query(Gene).filter_by(id=gid).first().ncbi_id for gid in annotatedGenes]

    ## remove genes covered that are not of the correct taxa (i.e. viral) 
    apQuery = [session.query(Uniprot).filter_by(id=uid).first() for uid in annotatedProts]
    _genesFromUniprot = [session.query(Gene).filter_by(id=uq.gene_id).first() for uq in apQuery]

    if not _genesFromUniprot:
        geneIdsFromUniprot = []
    else:
        remove_empty(_genesFromUniprot)
        genesFromUniprot = []

        for gene in _genesFromUniprot:
            if gene == None:
                continue
            geneTaxonId = session.query(Taxon).filter_by(id=gene.taxa_id).first()
            if int(geneTaxonId.ncbi_id) == int(taxonId):
                genesFromUniprot.append(gene)

        geneIdsFromUniprot = [str(g.id) for g in genesFromUniprot]

    annotatedGenes = list(set(annotatedGenes).union(set(geneIdsFromUniprot)))
    remove_empty(annotatedGenes)

    return [session.query(Gene).filter_by(id=gid).first().ncbi_id for gid in annotatedGenes]

def fetch_annotations(identifiers,engine,aspect='biological_process',
                      idType='uniprot',useIea=True,verbose=False):
    """
    Fetch the go annotations for a given list of identifiers

    If the identifier is 'uniprot'. Then combine the annotations from that uniportId 
    and the associated geneId if persent.

    If the identifier is 'geneid' the find all uniprot entries and combine all 
    uniprot and geneId results.

    The arg 'asTerms' return

    aspect is 'biological_process', 'cellular_component' or 'molecular_function'

    """

    acceptedCodes = get_evidence_codes(useIea=useIea)
    conn = engine.connect()

    if aspect not in ['biological_process','cellular_component','molecular_function']:
        raise Exception("Invalid aspect specified")

    ## error check
    if type(identifiers) != type([]):
        raise Exception("Takes a list of identifiers")

    annotations = {}
    taxaList = set([])
    idType = idType.lower()
    if idType not in ['uniprot','ncbi']:
        raise Exception("Invalid idType argument in fetch annotations use 'uniprot' or 'ncbi'")

    if idType == 'ncbi':
        if verbose:
            print('fetching annotations')
            print('...translating gene queries')
        s = select([Gene.id,Gene.ncbi_id]).where(Gene.ncbi_id.in_(identifiers))
        _geneQueries = conn.execute(s)
        geneQueries = _geneQueries.fetchall()
        if verbose:
            print("...%s/%s gene queries present"%(len(geneQueries),len(list(set(identifiers)))))
        
        for gquery in geneQueries:
            ncbiId = str(gquery['ncbi_id'])
            annotations[ncbiId] = set([])
            
            ## add results from the gene id
            s = select([GoTerm.go_id,GoTerm.name],GoAnnotation.go_term_id==GoTerm.id).\
                where(GoAnnotation.gene_id == gquery['id']).\
                where(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                where(GoTerm.aspect==aspect)

            _results = conn.execute(s)
            results = [tuple(map(str,r)) for r in _results.fetchall()]
            
            if len(results) > 1:
                annotations[ncbiId].update(results)

            ## check for corresponding uniprot ids
            s = select([Uniprot.id,Uniprot.uniprot_ac,Uniprot.gene_id]).\
                where(Uniprot.gene_id==gquery['id'])
            _uniprotQueries = conn.execute(s)
            uniprotQueries = _uniprotQueries.fetchall()
       
            if len(uniprotQueries) > 0:
                for uquery in uniprotQueries:
                    uniprotAc = str(uquery['uniprot_ac'])
            
                ## add results from the uniprot id
                s = select([GoTerm.go_id,GoTerm.name],GoAnnotation.go_term_id==GoTerm.id).\
                    where(GoAnnotation.uniprot_id == uquery['id']).\
                    where(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                    where(GoTerm.aspect==aspect)
              
                _results = conn.execute(s)
                results = [tuple(map(str,r)) for r in _results.fetchall()]
            
                if len(results) > 1:
                    annotations[ncbiId].update(results)
    
    elif idType == 'uniprot':
        if verbose:
            print('fetching annotations')
            print('...translating uniprot queries')
        s = select([Uniprot.id,Uniprot.uniprot_ac,Uniprot.gene_id]).where(Uniprot.uniprot_ac.in_(identifiers))
        _uniprotQueries = conn.execute(s)
        uniprotQueries = _uniprotQueries.fetchall()
        if verbose:
            print("...%s/%s uniprot queries present"%(len(uniprotQueries),len(list(set(identifiers)))))

        for uquery in uniprotQueries:
            uniprotAc = str(uquery['uniprot_ac'])
            annotations[uniprotAc] = set([])
            
            ## add results from the uniprot id
            s = select([GoTerm.go_id,GoTerm.name],GoAnnotation.go_term_id==GoTerm.id).where(GoAnnotation.uniprot_id == uquery['id']).\
                where(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                where(GoTerm.aspect==aspect)
              
            _results = conn.execute(s)
            results = [tuple(map(str,r)) for r in _results.fetchall()]
            
            if len(results) > 1:
                annotations[uniprotAc].update(results)

            ## add results from the associated gene id
            if uquery['gene_id']:
                s = select([GoTerm.go_id,GoTerm.name],GoAnnotation.go_term_id==GoTerm.id).\
                    where(GoAnnotation.gene_id == uquery['gene_id']).\
                    where(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                    where(GoTerm.aspect==aspect)
              
                _results = conn.execute(s)
                results = [tuple(map(str,r)) for r in _results.fetchall()]

                if len(results) > 1:
                    annotations[uniprotAc].update(results)
        
    ## remove any null results
    for key,items in annotations.iteritems():
        while () in items:
            items.remove(())

        annotations[key] = list(items)

    return annotations


def fetch_taxa_annotations(identifiers,engine,aspect='biological_process',
                           useIea=True,verbose=False):
    """
    Fetch the go annotations for a given list of taxa

    If the identifier is 'uniprot'. Then combine the annotations from that uniportId 
    and the associated geneId if persent.

    If the identifier is 'geneid' the find all uniprot entries and combine all 
    uniprot and geneId results.

    aspect is 'biological_process', 'cellular_component' or 'molecular_function'

    """

    acceptedCodes = get_evidence_codes(useIea=useIea)
    conn = engine.connect()

    if aspect not in ['biological_process','cellular_component','molecular_function']:
        raise Exception("Invalid aspect specified")

    ## error check
    if type(identifiers) != type([]):
        raise Exception("Takes a list of identifiers")

    goTerms = {}
    geneAnnotations = {}
    uniprotAnnotations = {}
    gene2uniprot,uniprot2gene = {},{}
    taxaList = set([])

    ## translate the ncbi taxa ids into db ids
    if verbose:
        print('fetching annotations')
        print('...translating uniprot queries')
    s = select([Taxon.id,Taxon.ncbi_id]).where(Taxon.ncbi_id.in_(identifiers))
    _taxaQueries = conn.execute(s)
    taxaQueries = _taxaQueries.fetchall()
    if verbose:
        print("...%s/%s taxa queries present"%(len(taxaQueries),len(list(set(identifiers)))))

    for tquery in taxaQueries:
        ## get all the annotations for the taxa
        s = select([GoTerm.go_id,GoTerm.name,GoAnnotation.uniprot_id,GoAnnotation.gene_id],GoAnnotation.taxa_id==tquery['id']).\
                    where(GoAnnotation.go_term_id == GoTerm.id).\
                    where(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                    where(GoTerm.aspect==aspect)
        _annotQueries = conn.execute(s)
        annotQueries = _annotQueries.fetchall()
    
        ## get all the gene ids and uniprot ids
        geneIds = []
        uniprotIds = []
        for aq in annotQueries:
            if aq['gene_id']:
                geneIds.append(aq['gene_id'])
            if aq['uniprot_id']:
                uniprotIds.append(aq['uniprot_id'])
            goTerms[str(aq[0])] = str(aq[1])

        ## get the gene to uniprot mappings
        s = select([Gene.id,Gene.ncbi_id],Gene.taxa_id==tquery['id'])
        _geneQueries = conn.execute(s)
        gene2Id = dict([tuple(map(str,r)) for r in _geneQueries.fetchall()])
        
        uniprot2id = {}
        if len(uniprotIds) > 0:
            s = select([Uniprot.id,Uniprot.uniprot_ac,Uniprot.gene_id]).where(Uniprot.id.in_(uniprotIds))
            _uniprotQueries = conn.execute(s)
            uniprotQueries = _uniprotQueries.fetchall()
            for uquery in uniprotQueries:
                uniprotAc = str(uquery['uniprot_ac'])
                uniprot2id[str(uquery['id'])] = uniprotAc
                if uquery['gene_id']:
                    if not gene2uniprot.has_key(str(uquery['gene_id'])):
                        if gene2Id.has_key(str(uquery['gene_id'])):
                            geneNcbi = gene2Id[str(uquery['gene_id'])] 
                            gene2uniprot[geneNcbi] = []
                        else:
                            geneNcbi = None
                    if geneNcbi:
                        gene2uniprot[geneNcbi].append(str(uquery['uniprot_ac']))
                        uniprot2gene[uniprotAc] = geneNcbi
                    else:
                        print "could not find", uquery['gene_id'], 'but we have', uniprotAc

        ## gene centric results
        for aq in annotQueries:
            result = [str(aq[0])]
            if aq[3]:
                geneNcbiId = gene2Id[str(aq[3])]
                if not geneAnnotations.has_key(geneNcbiId):
                    geneAnnotations[geneNcbiId] = set([])
                geneAnnotations[geneNcbiId].update(result)
            if aq[2]:
                uniprotAc = uniprot2id[str(aq[2])]
                if not uniprotAnnotations.has_key(uniprotAc):
                    uniprotAnnotations[uniprotAc] = set([])
                uniprotAnnotations[uniprotAc].update(result)
                
            ## check if we can map gene annotation to uniprot products
            if aq[3]:
                if gene2uniprot.has_key(geneNcbiId) and len(gene2uniprot[geneNcbiId]) > 0:
                    for newUniprotAc in gene2uniprot[geneNcbiId]:
                        if not uniprotAnnotations.has_key(newUniprotAc):
                            uniprotAnnotations[newUniprotAc] = set([])
                        uniprotAnnotations[newUniprotAc].update(result)
            
            ## check if we can map uniprot annotation to gene products
            if aq[2]:
                if uniprot2gene.has_key(uniprotAc):
                    newGeneId = uniprot2gene[uniprotAc]
                    if not geneAnnotations.has_key(newGeneId):
                        geneAnnotations[newGeneId] = set([])
                    geneAnnotations[newGeneId].update(result)

    ## prep the results
    for annotations in [geneAnnotations,uniprotAnnotations]:
        for key,items in annotations.iteritems():
            annotations[key] = list(items)

    return geneAnnotations,uniprotAnnotations

def read_annotation_file():
    """
    read the annotation file into a dictionary
    This will take some time
    This function is intended for use with database population

    http://www.geneontology.org/GO.format.gaf-2_0.shtml
    """

    annotationFile = get_annotation_file()
    annotationFid = open(annotationFile,'rU')
    result = {}

    for record in annotationFid:
        record = record[:-1].split("\t")

        ## check that it is a uniprot entry
        if record[0][0] == "!":
            continue
        
        dbObjectId = record[1]
        dbObjectSymbol = record[2]
        goID = record[4]
        dbReference = record[5]
        evidenceCode = record[6]
        aspect = record[8]
        dbObjectType = record[11]
        taxon = re.sub("taxon:","",record[12])
        date = record[13]
        assignedBy = record[14]

        ## ignore annotations with multiple species
        if re.search("\|",taxon):
            continue

        if not result.has_key(dbObjectId):
            result[dbObjectId] = {'names':set([]),'annots':{},'taxon':taxon}

        result[dbObjectId]['annots'][goID] = [aspect,evidenceCode]
        result[dbObjectId]['names'].update([dbObjectSymbol])

    annotationFid.close()
    return result
