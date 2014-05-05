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

def read_gene_info_file(geneIds=None,short=False):
    """
    read the essential info from NCBI's gene info file
    """

    geneInfoFile = os.path.join(CONFIG['data'],"gene_info.db")
    geneInfoFid = open(geneInfoFile,'rU')
    taxaList = set([])
    header = geneInfoFid.next()
    geneInfo ={}

    for record in geneInfoFid:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]):
            continue

        ncbiId = record[1]
        if geneIds == None:
            pass
        elif geneIds.has_key(ncbiId) == False:
            continue

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
    return geneInfo

def get_geneids_from_idmapping():
    """
    return a unique dict of geneids with refseq values from the idmapping file
    """

    geneInfo = read_gene_info_file(short=True)
    print 'loaded gene info in get_geneids...'

    ## get all unique geneIds in idmapping file
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rU')
    totalLines = 0
    lineCount = 0
    geneIds = {}

    for record in idmappingFid:
        record = record[:-1].split("\t")
        ncbiId = record[2]
        lineCount += 1
        
        if ncbiId == '':
            continue

        if not geneInfo.has_key(ncbiId):
            _geneIds = [re.sub("\s+","",_ncid) for _ncid in ncbiId.split(";")]
            ncbiId = None
        
            for _gid in _geneIds:
                if geneInfo.has_key(_gid):
                    ncbiId = _gid

        if geneIds.has_key(ncbiId) or ncbiId == None:
            continue

        geneIds[ncbiId] = record[3]

    idmappingFid.close()
    return geneIds,lineCount

def populate_taxon_table(taxonList,engine):
    """
    given a list of taxon ids populate the taxon table    
    """
    
    total = len(taxonList)
    wayPoints = [round(int(w)) for w in np.linspace(0,total,100)]
    taxonList = list(set([str(tax) for tax in taxonList]))
    namesFile = os.path.join(CONFIG['data'],"names.dmp")
    if os.path.exists(namesFile) == False:
        raise Exception("Cannot find names.dmp... exiting")

    namesFID = open(namesFile,'rU')
    taxaCount = 0
    timeStart = time.time()
    toAdd = {}

    for linja in namesFID:
        linja = linja.rstrip("\n")
        linja = linja.split("|")
        linja = [re.sub("\t","",element) for element in linja]

        scientificName,commonName = None,None
        if linja[3] == 'scientific name':
            taxID = linja[0]
            scientificName = linja[1]
        elif re.search("common name",linja[3]):
            taxID = linja[0]
            commonName = linja[1]
        else:
            continue

        ## only populate a subset of the taxa
        if taxID not in taxonList:
            continue
        
        ## if record does not exist
        if not toAdd.has_key(taxID) and scientificName != None:
            taxaCount += 1
            someTaxon = Taxon(taxID,name=scientificName)
            toAdd[taxID] = {'ncbi_id':taxID,'name':scientificName,'common_name_1':'',
                            'common_name_2':'','common_name_3':''}
            
        ## if record exists add a common name 
        elif toAdd.has_key(taxID) and commonName != None:
            if  toAdd[taxID]['common_name_1'] == '':
                toAdd[taxID]['common_name_1'] = commonName
            elif  toAdd[taxID]['common_name_2'] == '':
                toAdd[taxID]['common_name_2'] = commonName
            elif  toAdd[taxID]['common_name_3'] == '':
                toAdd[taxID]['common_name_3'] = commonName
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

def populate_gene_table(geneIds,taxaList,engine):
    """
    use the geneids derived from the idmapping file along with gene_info data to populate the gene table 
    """

    geneInfo = read_gene_info_file(geneIds=geneIds)
    timeStart = time.time()
    toAdd = []
    totalRecords = 0
    total = len(geneIds)
    wayPoints = [round(int(w)) for w in np.linspace(0,total,10)]

    ## batch query for taxa of each gene id
    with engine.begin() as connection:
        taxaQueries = connection.execute(Taxon.__table__.select(Taxon.ncbi_id.in_(taxaList)))
    taxaIds = {}
    for tq in taxaQueries:
        taxaIds[str(tq.ncbi_id)] = tq.id

    for ncbiId,refseq in geneIds.iteritems():
        taxId,symbol,synonyms,description = geneInfo[ncbiId]

        ## determine if record exists and add common names up until 3
        taxa_id = taxaIds[taxId]

        ## define the table entry
        toAdd.append({'ncbi_id':ncbiId,'description':description,
                      'symbol':symbol,'synonyms':synonyms,'taxa_id':taxa_id})
        totalRecords += 1
        
        if len(toAdd) >= 10000:
            with engine.begin() as connection:
                connection.execute(Gene.__table__.insert().
                                   values(toAdd))
            toAdd = []

        ## show progress
        if totalRecords in wayPoints:
            print("\t%s / %s"%(totalRecords,total))

    print('committing changes...')            
    with engine.begin() as connection:
        connection.execute(Gene.__table__.insert().
                           values(toAdd))

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique genes were added."%totalRecords
    return timeStr,addedStr
    
def populate_uniprot_table(lineCount,session,engine):
    """
    populate the uniprot table with entries from idmappings
    """
    
    print("...getting gene info")
    geneInfo = read_gene_info_file(short=True)
    timeStart = time.time()
    toAdd = []
    totalRecords = 0
    wayPoints = [round(int(w)) for w in np.linspace(0,lineCount,20)]
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rU')

    print("...getting keys from Gene table")
    geneIdMap = {}
    for g in session.query(Gene).yield_per(5):
        geneIdMap[g.ncbi_id] = g.id
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
        elif not geneInfo.has_key(ncbiId):
            _geneIds = [re.sub("\s+","",_ncid) for _ncid in ncbiId.split(";")]
            ncbiId = None
        
            for _gid in _geneIds:
                if geneInfo.has_key(_gid):
                    ncbiId = _gid
        
        if ncbiId != '' and ncbiId != None:
            gene_id = str(ncbiId)
        else:
            gene_id = None
        
        toAdd.append({'uniprot_id':uniprotKbAc,'uniprot_entry':uniprotKbEntry,
                      'refseq':refseq,'uniprot_taxa_id':uniprotTaxon,'gene_id':gene_id})
        
        totalRecords += 1
    
        if len(toAdd) >= 1000:
            for ta in toAdd:
                if ta['gene_id'] == None:
                    continue
                if not geneIdMap.has_key(ta['gene_id']):
                    raise Exception("Gene %s not present -- try repopulating the Gene table"%ta['gene_id'])
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
        if ta['gene_id'] == None:
            continue
        if not geneIdMap.has_key(ta['gene_id']):
            raise Exception("Gene %s not present -- try repopulating the Gene table"%ta['gene_id'])
            ta['gene_id'] = geneIdMap[ta['gene_id']]

    with engine.begin() as connection:
        connection.execute(Uniprot.__table__.insert().
                           values(toAdd))
           
    idmappingFid.close()
    del geneIdMap
    del geneInfo
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
            goName,goNamespace,goDef = None,None,None
            isObsolete = False
            termCount += 1
            toAdd[goId] = {'go_id': goId,'aspect':None,'name':None,
                           'description':None}
            continue

        ## find namespace and description  
        if re.search("^name\:",linja):
            goName = re.sub("^name\:\s+","",linja)
        if re.search("^namespace\:",linja):
            goNamespace = re.sub("^namespace\:\s+","",linja)
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
    
    print('committing changes...')
    with engine.begin() as connection:
        connection.execute(GoTerm.__table__.insert().
                           values(toAdd.values()))

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique uniprot entries were added."%termCount
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
    result = {}

    print("...getting keys from Uniprot table")
    uniprotIdMap = {}
    for u in session.query(Uniprot).yield_per(5):
        uniprotIdMap[str(u.uniprot_id)] = u.id

    print("...getting keys from GoTerm table")
    termIdMap = {}
    for g in session.query(GoTerm).yield_per(5):
        termIdMap[g.go_id] = g.id

    print("...getting keys from Taxon table")
    taxaIdMap = {}
    for t in session.query(Taxon).yield_per(5):
        taxaIdMap[str(t.ncbi_id)] = t.id
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

        if len(toAdd) >= 10000:

            for ta in toAdd:
                ta['go_term_id'] = termIdMap[ta['go_term_id']]
                ta['uniprot_id'] = uniprotIdMap[ta['uniprot_id']]
                ta['taxa_id'] = taxaIdMap[ta['taxa_id']]

            with engine.begin() as connection:
                connection.execute(GoAnnotation.__table__.insert().
                                   values(toAdd))
            toAdd = []

    print('committing final changes...')
    for ta in toAdd:
        ta['go_term_id'] = termIdMap[ta['go_term_id']]
        ta['uniprot_id'] = uniprotIdMap[ta['uniprot_id']]
        ta['taxa_id'] = taxaIdMap[ta['taxa_id']]

    with engine.begin() as connection:
        connection.execute(GoAnnotation.__table__.insert().
                           values(toAdd))
    ## clean up
    annotationFid.close()
    del uniprotIdMap

    print("...getting keys from Gene table")
    geneIdMap = {}
    for g in session.query(Gene).yield_per(5):
        geneIdMap[str(g.gene_id)] = g.id

    ## add annotations from gene2go
    print("...getting annotations from gene2go")
    header = gene2goFid.next()
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
    
        toAdd.append({'go_term_id':goId,'evidence_code':evidenceCode,
                      'pubmed_refs':pubmedRefs,'uniprot_id':None,
                      'gene_id':ncbiId,'taxa_id':taxId})

        if len(toAdd) >= 10000:

            for ta in toAdd:
                ta['go_term_id'] = termIdMap[ta['go_term_id']]
                ta['gene_id'] = geneIdMap[ta['gene_id']]
                ta['taxa_id'] = taxaIdMap[ta['taxa_id']]

            with engine.begin() as connection:
                connection.execute(GoAnnotation.__table__.insert().
                                   values(toAdd))
            toAdd = []

    print('committing final changes...')
    for ta in toAdd:
        ta['go_term_id'] = termIdMap[ta['go_term_id']]
        ta['gene_id'] = geneIdMap[ta['gene_id']]
        ta['taxa_id'] = taxaIdMap[ta['taxa_id']]

    with engine.begin() as connection:
        connection.execute(GoAnnotation.__table__.insert().
                           values(toAdd))
    del taxaIdMap
    del termIdMap
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique uniprot entries were added."%annotationCount
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

def get_taxa_list():
    """
    using the annotation and the gene info file return a unique list of taxa ids
    """

    annotationFile = get_annotation_file()
    annotationFid = open(annotationFile,'rU')
    annotsCount = 0
    annotatedIds = {}
    totalAnnotations = 0
    taxaList = set([])

    for record in annotationFid:
        record = record[:-1].split("\t")
        if record[0][0] == "!":
            continue
        if record[0] != 'UniProtKB':
            continue

        annotsCount += 1
        taxon = re.sub("taxon:","",record[12])
        if taxon == "" or re.search("\|",taxon):
            continue

        annotatedIds[record[1]] = None
        taxaList.update([taxon])
        totalAnnotations += 1

    geneInfoFile = os.path.join(CONFIG['data'],"gene_info.db")
    geneInfoFid = open(geneInfoFile,'rU')
    header = geneInfoFid.next()
    for record in geneInfoFid:
        record = record.rstrip("\n")
        record = record.split("\t")
        if re.search("^\#",record[0]):
            continue
        taxaList.update([record[0]])

    annotationFid.close()
    geneInfoFid.close()

    return list(taxaList),totalAnnotations

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
