#!/usr/bin/env python
"""
These are helper scripts to populate the database
"""

### make imports
import sys,os,re,time,csv
import sqlalchemy
import getpass
import numpy as np
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from htsint.database import get_annotation_file, get_ontology_file, get_gene2go_file
from DatabaseTables import Base,Taxon,Gene,Uniprot,GoTerm,GoAnnotation
from DatabaseTables import taxa_mapper,gene_mapper,uniprot_mapper,goterm_mapper

from htsint import __basedir__
sys.path.append(__basedir__)
try:
    from configure import CONFIG
except:
    CONFIG = None

def check_version():
    """
    ensure that sqlalchemy is at least 0.8
    """
    if float(sqlalchemy.__version__[:-2]) < 0.8:
        print "ERROR: SQLAlchemy version is less than required version"
        sys.exit()

def ask_upass():
    """
    returns the pass word for the database
    """

    if CONFIG == None:
        raise Exception("You must create a configure.py file and database before using database functions")
    
    check_version()
    upass = CONFIG['dbpass']
    if upass == '':
        upass = getpass.getpass()

    return upass

def db_connect(verbose=False,upass=''):
    """
    generic function to connect to db

    """

    if CONFIG == None:
        raise Exception("You must create a configure.py file and database before using database functions")
    
    check_version()

    ## declare variables
    uname = CONFIG['dbuser']
    dbhost = CONFIG['dbhost']
    dbname = CONFIG['dbname']
    port = CONFIG['dbport']

    ## get data base parameters
    if upass == '':
        upass = ask_upass()

    if dbname == '' or port == '' or dbhost == '' or uname=='':
        raise Exception("Invalid database parameters")
            
    ## create connection to db and create necessary tables 
    print "connecting to database: %s"%dbname
    engine = create_engine('postgres://%s:%s@%s:%s/%s'%(uname,upass,dbhost,port,dbname),echo=verbose)
    connection = engine.connect()
    Session = sessionmaker(bind=engine)
    session = Session()
    print 'connected.'

    return session,engine

def read_gene_info_file(lineCount=False,short=False):
    """
    read the essential info from NCBI's gene info file
    """

    geneInfoFile = os.path.join(CONFIG['data'],"gene_info.db")
    geneInfoFid = open(geneInfoFile,'rU')
    header = geneInfoFid.next()
    geneInfo ={}
    totalLines = 0

    for record in geneInfoFid:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]):
            continue

        if lineCount == True:
            totalLines += 1
            continue

        ncbiId = record[1]
        taxId = record[0]
        symbol = record[2]
        synonyms = record[4]
        chromosome = record[6]
        map_location = record[7]
        description = record[8]
        if short == False:
            geneInfo[ncbiId] = [taxId,symbol,synonyms,description]
        else:
            geneInfo[ncbiId] = taxId

    geneInfoFid.close()
    
    if lineCount == True:
        return totalLines
    else:
        return geneInfo

def get_file_sizes():
    """
    return a unique dict of geneids with refseq values from the idmapping file
    """

    geneInfoCount = read_gene_info_file(lineCount=True)
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rU')
    lineCount = 0

    for record in idmappingFid:
        lineCount += 1
        
    idmappingFid.close()
    return lineCount,geneInfoCount

def populate_taxon_table(engine):
    """
    given a list of taxon ids populate the taxon table    
    """
    
    namesFile = os.path.join(CONFIG['data'],"names.dmp")
    if os.path.exists(namesFile) == False:
        raise Exception("Cannot find names.dmp... exiting")

    namesFID = open(namesFile,'rU')
    taxaCount = 0
    timeStart = time.time()
    toAdd = {}
    taxaID = None
    debug = 0

    for linja in namesFID:
        debug += 1

        linja = linja.rstrip("\n")
        linja = linja.split("|")
        linja = [re.sub("\t","",element) for element in linja]

        scientificName,commonName = None,None
        if linja[3] == 'scientific name':
            taxaID = linja[0]
            scientificName = linja[1]
        elif re.search("common name",linja[3]):
            taxaID = linja[0]
            commonName = linja[1]
        else:
            continue

        if taxaID in ['root']:
            continue
                
        ## if record does not exist
        if not toAdd.has_key(taxaID):
            taxaCount += 1
            if len(toAdd) >= 300000:
                with engine.begin() as connection:
                    connection.execute(Taxon.__table__.insert().
                                       values(toAdd.values()))
                toAdd = {}

            toAdd[taxaID] = {'ncbi_id':taxaID,'name':None,'common_name_1':None,
                            'common_name_2':None,'common_name_3':None}
        
        ## if record exists add a common name 
        if toAdd.has_key(taxaID) and scientificName != None:
            toAdd[taxaID]['name'] = scientificName
        elif toAdd.has_key(taxaID) and commonName != None:
            if  toAdd[taxaID]['common_name_1'] == None:
                toAdd[taxaID]['common_name_1'] = commonName
            elif  toAdd[taxaID]['common_name_2'] == None:
                toAdd[taxaID]['common_name_2'] = commonName
            elif  toAdd[taxaID]['common_name_3'] == None:
                toAdd[taxaID]['common_name_3'] = commonName
        else:
            continue

    print('committing changes...')
    with engine.begin() as connection:
        connection.execute(Taxon.__table__.insert().
                           values(toAdd.values()))
    namesFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr =  "...%s unique taxa were added."%taxaCount
    return timeStr, addedStr

def populate_gene_table(geneInfoCount,session,engine):
    """
    use the geneids derived from the idmapping file along with gene_info data to populate the gene table 
    """

    timeStart = time.time()
    toAdd = []
    totalRecords = 0
    total = geneInfoCount
    wayPoints = [round(int(w)) for w in np.linspace(0,total,20)]
    geneInfoFile = os.path.join(CONFIG['data'],"gene_info.db")
    geneInfoFid = open(geneInfoFile,'rU')
    header = geneInfoFid.next()
    taxaIdMap = taxa_mapper(session)

    for record in geneInfoFid:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]):
            continue

        taxId = record[0]
        ncbiId = record[1]
        symbol = record[2]
        synonyms = record[4]
        #chromosome = record[6]
        #map_location = record[7]
        description = record[8]
       
        ## define the table entry
        toAdd.append({'ncbi_id':ncbiId,'description':description,
                      'symbol':symbol,'synonyms':synonyms,'taxa_id':taxId})
        totalRecords += 1
        
        if len(toAdd) >= 300000:
            toRemove = []
            for ta in toAdd:
                if taxaIdMap.has_key(ta['taxa_id']):
                    ta['taxa_id'] = taxaIdMap[ta['taxa_id']]
                else:
                    toRemove.append(ta)

            for ta in toRemove:
                toAdd.remove(ta)

            with engine.begin() as connection:
                connection.execute(Gene.__table__.insert().
                                   values(toAdd))
            toAdd = []

        ## show progress
        if totalRecords in wayPoints:
            print("\t%s / %s"%(totalRecords,total))

    print('committing changes...')
    toRemove = []
    for ta in toAdd:
        if taxaIdMap.has_key(ta['taxa_id']):
            ta['taxa_id'] = taxaIdMap[ta['taxa_id']]
        else:
            toRemove.append(ta)

    for ta in toRemove:
        toAdd.remove(ta)

    del taxaIdMap

    with engine.begin() as connection:
        connection.execute(Gene.__table__.insert().
                           values(toAdd))

    geneInfoFid.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique genes were added."%totalRecords
    return timeStr,addedStr
    
def populate_uniprot_table(lineCount,session,engine):
    """
    populate the uniprot table with entries from idmappings
    """
    
    print("...getting gene info")
    geneIdMap = gene_mapper(session)
    timeStart = time.time()
    toAdd = []
    totalRecords = 0
    wayPoints = [round(int(w)) for w in np.linspace(0,lineCount,20)]
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rU')
    print("...populating rows")

    for record in idmappingFid:
        record = record[:-1].split("\t")

        uniprotKbAc = record[0]
        uniprotKbEntry = record[1]
        ncbiId = record[2]
        refseq = record[3]
        uniprotTaxon = record[13]

        if ncbiId == '':
            pass
        elif not geneIdMap.has_key(ncbiId):
            _geneIds = [re.sub("\s+","",_ncid) for _ncid in ncbiId.split(";")]
            ncbiId = None
        
            for _gid in _geneIds:
                if geneIdMap.has_key(_gid):
                    ncbiId = _gid
        
        if ncbiId != '' and ncbiId != None:
            gene_id = str(ncbiId)
        else:
            gene_id = None
        
        toAdd.append({'uniprot_id':uniprotKbAc,'uniprot_entry':uniprotKbEntry,
                      'refseq':refseq,'uniprot_taxa_id':uniprotTaxon,'gene_id':gene_id})
        
        totalRecords += 1
    
        if len(toAdd) >= 100000:
            for ta in toAdd:
                if not geneIdMap.has_key(ta['gene_id']):
                    ta['gene_id'] = None
                if ta['gene_id'] == None:
                    continue
                
                ta['gene_id'] = geneIdMap[ta['gene_id']]

            with engine.begin() as connection:
                connection.execute(Uniprot.__table__.insert().
                                   values(toAdd))
            toAdd = []
    

        ## show progress
        if totalRecords in wayPoints:
            print("\t%s / %s"%(totalRecords,lineCount))

    print('committing final changes...')
    for ta in toAdd:
        if not geneIdMap.has_key(ta['gene_id']):
            ta['gene_id'] = None
        if ta['gene_id'] == None:
            continue
        
        ta['gene_id'] = geneIdMap[ta['gene_id']]

    del geneIdMap

    with engine.begin() as connection:
        connection.execute(Uniprot.__table__.insert().
                           values(toAdd))
           
    idmappingFid.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique uniprot entries were added."%totalRecords
    return timeStr,addedStr

def populate_go_terms(engine):
    """ 
    read in the obo file and use it to populate the terms
    """

    timeStart = time.time()
    ontologyFile = get_ontology_file()
    fid = open(ontologyFile,'r')
    termCount = 0
    goId,goName,goNamespace,goDef = None,None,None,None
    toAdd = {}
    isObsolete = False

    for linja in fid.readlines():
        linja = linja[:-1]

        ## find go id
        if re.search("^id\:",linja):
            goId = re.sub("^id\:|\s+","",linja)
            goName,goNamespace,goDef,goAltId = None,None,None,None
            isObsolete = False
            termCount += 1
            toAdd[goId] = {'go_id': goId,'aspect':None,'name':None,
                           'alternate_id':None,'description':None}
            continue

        ## find namespace and description  
        if re.search("^name\:",linja):
            goName = re.sub("^name\:\s+","",linja)
        if re.search("^namespace\:",linja):
            goNamespace = re.sub("^namespace\:\s+","",linja)
        if re.search("^alt_id\:",linja):
            goAltId = re.sub("^alt_id\:\s+","",linja)
        if re.search("^def\:",linja):
            goDef = re.sub("^def\:\s+","",linja)
            if re.search("OBSOLETE\.",goDef):
                isObsolete = True
        
        ## ignore obolete terms
        if isObsolete == True:
            goId = None

        if goId != None:
            if goNamespace != None and toAdd[goId]['aspect'] == None:
                toAdd[goId]['aspect'] = goNamespace
            if goDef != None and toAdd[goId]['description'] == None:
                toAdd[goId]['description'] = goDef
            if goName != None and toAdd[goId]['name'] == None:
                toAdd[goId]['name'] = goName
            if goAltId != None and toAdd[goId]['alternate_id'] == None:
                toAdd[goId]['alternate_id'] = goAltId

    print('committing changes...')
    with engine.begin() as connection:
        connection.execute(GoTerm.__table__.insert().
                           values(toAdd.values()))

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique go term entries were added."%termCount
    return timeStr,addedStr

def populate_go_annotations(totalAnnotations,session,engine):
    """
    read the annotation file into a dictionary
    This will take some time
    This function is intended for use with                                                    
    http://www.geneontology.org/GO.format.gaf-2_0.shtml
    """

    timeStart = time.time()
    toAdd = []
    annotationFile = get_annotation_file()
    annotationFid = open(annotationFile,'rU')
    gene2goFile = get_gene2go_file()
    gene2goFid = open(gene2goFile,'rU')
    wayPoints = [round(int(w)) for w in np.linspace(0,totalAnnotations,20)]
    annotationCount = 0

    print("...loading mappers")
    termIdMap = goterm_mapper(session)
    taxaIdMap = taxa_mapper(session)
    uniprotIdMap = uniprot_mapper(session)
    print("...populating rows")

    ## add annotations from uniprot annotation file
    print("...getting annotations from gene_association (uniprot)")
    for record in annotationFid:
        record = record[:-1].split("\t")

        ## check that it is a uniprot entry
        if record[0][0] == "!":
            continue
        if record[0] != 'UniProtKB':
            continue
        
        dbObjectId = record[1]
        dbObjectSymbol = record[2]
        goId = record[4]
        pubmedRefs = record[5]
        evidenceCode = record[6]
        aspect = record[8]
        goTermName = record[11]
        taxon = re.sub("taxon:","",record[12])
        date = record[13]
        assignedBy = record[14]

        ## ignore annotations with multiple species
        if re.search("\|",taxon):
            continue

        ## update progress
        annotationCount += 1
        if annotationCount in wayPoints:
            print("\t%s / %s"%(annotationCount,totalAnnotations))

        toAdd.append({'go_term_id':goId,'evidence_code':evidenceCode,
                      'pubmed_refs':pubmedRefs,'uniprot_id':dbObjectId,
                      'gene_id':None,'taxa_id':taxon})

        if len(toAdd) >= 100000:
            
            toRemove = []
            for ta in toAdd:
                ## remove invalid term ids
                if not termIdMap.has_key(ta['go_term_id']):
                    queryTerm = session.query(GoTerm).filter_by(alternate_id=ta['go_term_id']).first()
                    if queryTerm == None:
                        toRemove.append(ta)
                        continue
                    ta['go_term_id'] = queryTerm.id
                else:
                    ta['go_term_id'] = termIdMap[ta['go_term_id']]

                ## remove invalid uniprot ids
                if not uniprotIdMap.has_key(ta['uniprot_id']):
                    toRemove.append(ta)
                    continue
                    
                ta['uniprot_id'] = uniprotIdMap[ta['uniprot_id']]
            
                if taxaIdMap.has_key(ta['taxa_id']):
                    ta['taxa_id'] = taxaIdMap[ta['taxa_id']]
                else:
                    ta['taxa_id'] = None

            if len(toRemove) > 0:
                print("removing...%s invalid annotations"%len(toRemove))
            for ta in toRemove:
                toAdd.remove(ta)

            with engine.begin() as connection:
                connection.execute(GoAnnotation.__table__.insert().
                                   values(toAdd))
            toAdd = []

    print('committing final changes...')
    toRemove = []
    for ta in toAdd:
        ## remove invalid term ids
        if not termIdMap.has_key(ta['go_term_id']):
            queryTerm = session.query(GoTerm).filter_by(alternate_id=ta['go_term_id']).first()
            if queryTerm == None:
                toRemove.append(ta)
                continue
            ta['go_term_id'] = queryTerm.id
        else:
            ta['go_term_id'] = termIdMap[ta['go_term_id']]

        ## remove invalid uniprot ids
        if not uniprotIdMap.has_key(ta['uniprot_id']):
            toRemove.append(ta)
            continue

        ta['uniprot_id'] = uniprotIdMap[ta['uniprot_id']]
            
        if taxaIdMap.has_key(ta['taxa_id']):
            ta['taxa_id'] = taxaIdMap[ta['taxa_id']]
        else:
            ta['taxa_id'] = None
 
    if len(toRemove) > 0:
        print("removing...%s invalid annotations"%len(toRemove))
    for ta in toRemove:
        toAdd.remove(ta)

    del uniprotIdMap
    annotationFid.close()

    with engine.begin() as connection:
        connection.execute(GoAnnotation.__table__.insert().
                           values(toAdd))

    ## add annotations from gene2go
    print("...getting annotations from gene2go")
    header = gene2goFid.next()
    geneIdMap = gene_mapper(session)
    toAdd = []

    for record in gene2goFid:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]) or len(record) != 8:
            continue
    
        taxId = record[0]
        ncbiId = record[1]
        goId = record[2]
        evidenceCode = record[3]
        qualifier = record[4]
        go_term_description = record[5]
        pubmedRefs = record[6]
        go_aspect = record[7]
        annotationCount += 1

        if annotationCount in wayPoints:
            print("\t%s / %s"%(annotationCount,totalAnnotations))

        toAdd.append({'go_term_id':goId,'evidence_code':evidenceCode,
                      'pubmed_refs':pubmedRefs,'uniprot_id':None,
                      'gene_id':ncbiId,'taxa_id':taxId})

        if len(toAdd) >= 100000:
            toRemove = []
            for ta in toAdd:
                ## remove invalid term ids
                if not termIdMap.has_key(ta['go_term_id']):
                    queryTerm = session.query(GoTerm).filter_by(alternate_id=ta['go_term_id']).first()
                    if queryTerm == None:
                        toRemove.append(ta)
                        continue
                    ta['go_term_id'] = queryTerm.id
                else:
                    ta['go_term_id'] = termIdMap[ta['go_term_id']]
            
                ## remove invalid gene ids
                if not geneIdMap.has_key(ta['gene_id']):
                    toRemove.append(ta)
                    continue

                ta['gene_id'] = geneIdMap[ta['gene_id']]

                if taxaIdMap.has_key(ta['taxa_id']):
                    ta['taxa_id'] = taxaIdMap[ta['taxa_id']]
                else:
                    ta['taxa_id'] = None

            if len(toRemove) > 0:
                print("removing...%s invalid annotations"%len(toRemove))
            for ta in toRemove:
                toAdd.remove(ta)

            with engine.begin() as connection:
                connection.execute(GoAnnotation.__table__.insert().
                                   values(toAdd))
            toAdd = []

    print('committing final changes...')
    toRemove = []
    for ta in toAdd:
        ## remove invalid term ids
        if not termIdMap.has_key(ta['go_term_id']):
            queryTerm = session.query(GoTerm).filter_by(alternate_id=ta['go_term_id']).first()
            if queryTerm == None:
                toRemove.append(ta)
                continue
            ta['go_term_id'] = queryTerm.id
        else:
            ta['go_term_id'] = termIdMap[ta['go_term_id']]

        ## remove invalid gene ids
        if not geneIdMap.has_key(ta['gene_id']):
            toRemove.append(ta)
            continue

        ta['gene_id'] = geneIdMap[ta['gene_id']]
        
        if taxaIdMap.has_key(ta['taxa_id']):
            ta['taxa_id'] = taxaIdMap[ta['taxa_id']]
        else:
            ta['taxa_id'] = None

    if len(toRemove) > 0:
        print("removing...%s invalid annotations"%len(toRemove))
    for ta in toRemove:
        toAdd.remove(ta)

    del taxaIdMap
    del geneIdMap
    del termIdMap

    with engine.begin() as connection:
        connection.execute(GoAnnotation.__table__.insert().
                           values(toAdd))

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique go annotation entries were added."%annotationCount
    return timeStr,addedStr

def print_db_summary():
    """
    print a summary of rows and tables for the database
    """
    
    session,engine = db_connect(verbose=False)
    print("\nDATABASE - %s - SUMMARY"%CONFIG['dbname'])
    for table in [Taxon,Gene,Uniprot,GoTerm,GoAnnotation]:
        print("There are %s entries in the %s table"%(session.query(table).count(),table.__tablename__))
    print "\n"


def get_idmapping_file():
    """
    check for presence of the annotation file
    raise exception when not found
    return the file path 
    """

    if CONFIG == None:
        raise Exception("You must create a configure.py before GeneOntology")

    dataDir = CONFIG['data']
    idmappingFile = os.path.join(dataDir,'idmapping.tb.db')
    if os.path.exists(idmappingFile) == False:
        raise Exception("Could not find 'idmapping.tb.db' -- did you run FetchDbData.py?")

    return idmappingFile

def print_go_summary(outfile=os.path.join(".","go_summary.csv")):
    """
    print a summary about the organisms with gene ontology information
    """

    session,engine = db_connect()
    goTaxa = []
    fid = open(outfile,'w')
    writer = csv.writer(fid)
    writer.writerow(["ncbi_id","name","common_name","total_genes","total_annotations"])

    for taxa in goTaxa:
        taxQuery = session.query(Taxon).filter_by(ncbi_id=taxa).first()
        geneQuery = session.query(Gene).filter_by(taxa_id=taxQuery.id)
        annotQuery = session.query(GoAnnotation).filter_by(taxa_id=taxQuery.id)
        
        print "..."
        print "TaxID:       %s"%taxQuery.ncbi_id
        print "Species:     %s"%taxQuery.name
        print "Common Name 1: %s"%taxQuery.common_name_1
        print "Common Name 2: %s"%taxQuery.common_name_2
        print "Common Name 3: %s"%taxQuery.common_name_3
        print "Num. Genes:  %s"%geneQuery.count()
        print "Num. GO Annotations:  %s"%annotQuery.count()

        if taxQuery.common_name_1 == '':
            taxQuery.common_name_1 = 'None'

        writer.writerow([taxQuery.ncbi_id,taxQuery.name,taxQuery.common_name_1,
                         geneQuery.count(),annotQuery.count()])

    fid.close()
