#!/usr/bin/env python
"""
library of functions for use with the GeneOntology class
"""

import os,sys,re
from htsint import __basedir__
from DatabaseTables import Taxon,Uniprot,Gene,GoTerm,GoAnnotation

sys.path.append(__basedir__)
try:
    from configure import CONFIG
except:
    CONFIG = None

def get_ontology_file():
    """
    check for presence of ontology file
    raise exception when not found
    return the file path
    """

    if CONFIG == None:
        raise Exception("You must create a configure.py before GeneOntology")

    dataDir = CONFIG['data']

    ontologyFile = os.path.join(dataDir,'go.obo')
    if os.path.exists(ontologyFile) == False:
        raise Exception("Could not find 'go.obo' -- did you run FetchDbData.py?")

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

    if CONFIG == None:
        raise Exception("You must create a configure.py before GeneOntology")

    dataDir = CONFIG['data']
    annotationFile = os.path.join(dataDir,fileName)
    if os.path.exists(annotationFile) == False:
        raise Exception("Could not find '%s' -- did you run FetchDbData.py?"%(fileName))

    return annotationFile

def get_total_annotations():
    """
    get the number of annotations in the uniprot file
    this does not include the gene2go file
    """

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

    if CONFIG == None:
        raise Exception("You must create a configure.py before GeneOntology")

    dataDir = CONFIG['data']
    annotationFile = os.path.join(dataDir,'gene2go.db')
    if os.path.exists(annotationFile) == False:
        raise Exception("Could not find 'gene2go' -- did you run FetchDbData.py?")

    return annotationFile

def get_evidence_codes():
    """
    returns the list of evidence codes to be used by default
    """

    expEvidCodes = ["EXP","IDA","IPI","IMP","IGI","IEP"]
    compEvidCodes = ["ISS","ISO","ISA","ISM","IGC","RCA"]
    statEvidCodes = ["TAS","NAS","IC"]
    nonCuratedEvidCodes = ["IEA"]

    return expEvidCodes + statEvidCodes + compEvidCodes

def fetch_annotations(identifiers,session,aspect='biological_process',
                      idType='uniprot',asTerms=True):
    """
    Fetch the go annotations for a given list of identifiers

    If the identifier is 'uniprot'. Then combine the annotations from that uniportId 
    and the associated geneId if persent.

    If the identifier is 'geneid' the find all uniprot entries and combine all 
    uniprot and geneId results.

    The arg 'asTerms' return

    aspect is 'biological_process', 'cellular_component' or 'molecular_function'

    """

    acceptedCodes = get_evidence_codes()

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
        geneQueries = session.query(Gene).filter(Gene.ncbi_id.in_(identifiers)).all()
                
        for geneQuery in geneQueries:
            taxaList.update([geneQuery.taxa_id])
            annotations[geneQuery.ncbi_id] = set([])
            
            ## add annotations from the ncbi gene id
            results = session.query(GoAnnotation).join(GoTerm).\
                      filter(GoAnnotation.gene_id==geneQuery.id).\
                      filter(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                      filter(GoTerm.aspect==aspect).all()
            if results:
                annotations[geneQuery.ncbi_id].update(results)

            ## add annotations from all associated uniprot ids
            uniprotQueries = session.query(Uniprot).filter_by(gene_id=geneQuery.id).all()
            for uniprotQuery in uniprotQueries:
                results = session.query(GoAnnotation).join(GoTerm).\
                          filter(GoAnnotation.uniprot_id==uniprotQuery.id).\
                          filter(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                          filter(GoTerm.aspect==aspect).all()
                if results:
                    annotations[geneQuery.ncbi_id].update(results)
    
    elif idType == 'uniprot':
        uniprotQueries = session.query(Uniprot).filter(Uniprot.uniprot_ac.in_(identifiers)).all()

        for uniprotQuery in uniprotQueries:
            annotations[uniprotQuery.uniprot_ac] = set([])
            
            ## add results from the uniprot id
            results = session.query(GoAnnotation).join(GoTerm).\
                      filter(GoAnnotation.uniprot_id==uniprotQuery.id).\
                      filter(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                      filter(GoTerm.aspect==aspect).all()
            if results:
                annotations[uniprotQuery.uniprot_ac].update(results)

            ## add results from the associated gene id
            geneQuery = session.query(Gene).filter_by(id=uniprotQuery.gene_id).first()
            if geneQuery != None:
                taxaList.update([geneQuery.taxa_id])
                results = session.query(GoAnnotation).join(GoTerm).\
                          filter(GoAnnotation.gene_id==geneQuery.id).\
                          filter(GoAnnotation.evidence_code.in_(acceptedCodes)).\
                          filter(GoTerm.aspect==aspect).all()
            if results:
                annotations[uniprotQuery.uniprot_ac].update(results)
    
    ## remove any null results
    for key,items in annotations.iteritems():
        annotations[key] = list(items)
        if None in items:
            annotations[key].remove(None)
        if asTerms == True:
            annotations[key] = [session.query(GoTerm).filter_by(id = a.go_term_id).first().go_id for a in items]

    return annotations

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
