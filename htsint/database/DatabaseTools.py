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
from htsint.database import get_annotation_file, get_ontology_file
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

def read_gene_info_file(geneIds):
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
        if geneIds.has_key(ncbiId) == False:
            continue

        taxId = record[0]
        symbol = record[2]
        synonyms = record[4]
        chromosome = record[6]
        map_location = record[7]
        description = record[8]
        geneInfo[ncbiId] = [taxId,symbol,synonyms,description]

    geneInfoFid.close()
    return geneInfo


def get_geneids_from_idmapping():
    """
    return a unique dict of geneids with refseq values from the idmapping file
    
    """

    ## get all unique geneIds in idmapping file
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rU')
    totalLines = 0
    lineCount = 0
    geneIds = {}
    for record in idmappingFid:
        record = record[:-1].split("\t")
        lineCount += 1
        if record[2] != '' and not geneIds.has_key(record[2]):
            geneIds[record[2]] = record[3]

    idmappingFid.close()
    return geneIds,lineCount

def populate_taxon_table(taxonList,session,engine):
    """
    given a list of taxon ids populate the taxon table
    
    """
    
    total = len(taxonList)
    wayPoints = [round(int(w)) for w in np.linspace(0,total,100)]
    taxonList = list(set([str(tax) for tax in taxonList]))
    namesFile = os.path.join(CONFIG['data'],"names.dmp")
    if os.path.exists(namesFile) == False:
        print "ERROR: Cannot find names.dmp... exiting"
        sys.exit()

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

            if len(toAdd) > 1000:
                with session.begin(subtransactions=True):
                    session.add_all(toAdd.values())
                toAdd = {}

            someTaxon = Taxon(taxID,name=scientificName)
            toAdd[taxID] = (someTaxon)
            
            ## show progress
            if taxaCount in wayPoints:
                print("\t%s percent finished"%(round(taxaCount/float(total)*100.0)))

        ## if record exists add a common name 
        elif toAdd.has_key(taxID) and commonName != None:
            if  toAdd[taxID].common_name_1 == '':
                toAdd[taxID].common_name_1 = commonName
            elif commonName not in [toAdd[taxID].common_name_1] and toAdd[taxID].common_name_2 == '':
                toAdd[taxID].common_name_2 = commonName
            elif commonName not in [toAdd[taxID].common_name_1,toAdd[taxID].common_name_2] and toAdd[taxID].common_name_3 == '':
                toAdd[taxID].common_name_3 = commonName
        else:
            continue

    print('committing changes...')
    with session.begin(subtransactions=True):
        session.add_all(toAdd.values())
    session.commit()
    namesFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr =  "...%s unique taxa were added."%taxaCount
    return timeStr, addedStr

def populate_gene_table(geneIds,session):
    """
    use the geneids derived from the idmapping file along with gene_info data to populate the gene table 

    """

    geneInfo = read_gene_info_file(geneIds)
    timeStart = time.time()
    toAdd = []
    totalRecords = 0
    total = len(geneIds)
    wayPoints = [round(int(w)) for w in np.linspace(0,total,100)]

    for ncbiId,refseq in geneIds.iteritems():
        
        if not geneInfo.has_key(ncbiId):
            "WARNING: key not found:  %s "%(ncbiId)
            continue

        taxId,symbol,synonyms,description = geneInfo[ncbiId]

        ## determine if record exists and add common names up until 3
        queryTax = session.query(Taxon).filter_by(ncbi_id=taxId).first()
        taxa_id = queryTax.id
        
        ## define the table entry
        someGene = Gene(ncbiId,description,symbol,synonyms,taxa_id)
        toAdd.append(someGene)
        totalRecords += 1
        
        ## show progress
        if totalRecords in wayPoints:
            print("\t%s percent finished"%(round(totalRecords/float(total))*100.0))

        ## periodically update the db
        if len(toAdd) > 1000:
            with session.begin(subtransactions=True):
                session.add_all(toAdd)
                toAdd = []

    ## add the remaining genes
    with session.begin(subtransactions=True):
        session.add_all(toAdd)
        toAdd = []
    session.commit()

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique genes were added."%totalRecords
    return timeStr,addedStr
    
def populate_uniprot_table(lineCount,session):
    """
    populate the uniprot table with entries from idmappings
    """
    
    timeStart = time.time()
    toAdd = []
    totalRecords = 0
    wayPoints = [round(int(w)) for w in np.linspace(0,lineCount,100)]
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rU')
    
    for record in idmappingFid:
        record = record[:-1].split("\t")

        uniprotKbAc = record[0]
        uniprotKbEntry = record[1]
        geneId = record[2]
        refseq = record[3]
        uniprotTaxon = record[13]
    
        queryGene = session.query(Gene).filter_by(ncbi_id=geneId).first()
        if queryGene == None:
            gene_id = None
        else:
            gene_id = queryGene.id

        uniprotEntry = Uniprot(uniprotKbAc,uniprotKbEntry,refseq,uniprotTaxon,gene_id)
        toAdd.append(uniprotEntry)
        totalRecords += 1

        ## show progress
        if totalRecords in wayPoints:
            print("\t%s percent finished"%(round(totalRecords/float(lineCount)*100.0)))

        ## periodically update the db
        if len(toAdd) > 1000:
            with session.begin(subtransactions=True):
                session.add_all(toAdd)
                session.flush()
                toAdd = []

    ## add the remaining entries
    session.add_all(toAdd)
    session.commit()

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique uniprot entries were added."%totalRecords
    return timeStr,addedStr


def populate_go_terms(session):
    """ 
    read in the obo file and use it to populate the terms
    """

    timeStart = time.time()
    ontologyFile = get_ontology_file()
    fid = open(ontologyFile,'r')
    termCount = 0
    goId = None
    toAdd = []
    isObsolete = False

    for linja in fid.readlines():
        linja = linja[:-1]

        ## find go id
        if re.search("^id\:",linja):
            goId = re.sub("^id\:|\s+","",linja)
            goName,goNamespace = None,None
            isObsolete = False
            termCount += 1
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

        ## add the term  
        if goId != None and re.search("^is\_a\:",linja):
            toAdd.append(GoTerm(goId,goNamespace,goDef))

        ## periodically update the db
        if len(toAdd) > 1000:
            with session.begin(subtransactions=True):
                session.add_all(toAdd)
                session.flush()
                toAdd = []

    ## add the remaining entries
    session.add_all(toAdd)
    session.commit()

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique uniprot entries were added."%termCount
    return timeStr,addedStr


def populate_go_annotations(session):
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
    annotationCount = 0
    result = {}

    for record in annotationFid:
        record = record[:-1].split("\t")

        ## check that it is a uniprot entry
        if record[0][0] == "!":
            continue
        if record[0] != 'UniProtKB':
            continue

        if not annotatedIds.has_key(record[1]):
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

        annotationCount += 1

        ## get the foreign keys
        queryTerm = session.query(GoTerm).filter_by(go_id=goId).first()
        go_term_id = queryTerm.id

        queryUniprot = session.query(Uniprot).filter_by(uniprot_id=dbObjectId).first()
        uniprot_id = queryUniprot.id

        queryTaxa = session.query(Taxon).filter_by(ncbi_id=taxon).first()
        taxa_id = queryTaxa.id

        toAdd.append(annotation(go_term_id,evidenceCode,pubmedRefs,uniprot_id,taxa_id))

        ## periodically update the db
        if len(toAdd) > 1000:
            with session.begin(subtransactions=True):
                session.add_all(toAdd)
                session.flush()
                toAdd = []

    ## add the remaining entries
    session.add_all(toAdd)
    session.commit()

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


## the following are deprecated


'''
def get_all_go_taxa():
    """
    get the list of all taxa that have gene ontology annotations
    """

    gene2goFile = os.path.join(CONFIG['data'],"gene2go.db")
    if os.path.exists(gene2goFile) == False:
        print "ERROR: populate_gene_table() exiting... could not find geneInfoFile"
        return

    gene2goFID = open(gene2goFile,'rU')
    header = gene2goFID.next()
    totalTerms,totalAnnotations = 0,0
    timeStart = time.time()
    toAddTerms = []
    toAddAnnotations = []
    allTaxa = set([])

    for record in gene2goFID:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]) or len(record) != 8:
            print 'continuing',len(record),record[0]
            continue

        allTaxa.update([record[0]])

    return list(allTaxa)
'''
