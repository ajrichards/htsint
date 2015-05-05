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
from htsint import Configure
from DatabaseTables import Base,Taxon,Gene,Uniprot,GoTerm,GoAnnotation
from DatabaseTables import taxa_mapper,gene_mapper,uniprot_mapper,goterm_mapper
from htsint.database import get_annotation_file, get_ontology_file, get_gene2go_file

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


    config = Configure()

    for key in ['data','dbname']:
        if config.log[key] == '':
            raise Exception("You must modify the config file before running DatabaseFetch.py")
    
    check_version()
    upass = config.log['dbpass']
    if upass == '':
        upass = getpass.getpass()

    return upass

def db_connect(verbose=False,upass=''):
    """
    generic function to connect to db

    """

    config = Configure()

    for key in ['data','dbname']:
        if config.log[key] == '':
            raise Exception("You must modify the config file before running DatabaseFetch.py")
    
    check_version()

    ## declare variables
    uname = config.log['dbuser']
    dbhost = config.log['dbhost']
    dbname = config.log['dbname']
    port = config.log['dbport']

    ## get data base parameters
    if upass == '':
        upass = ask_upass()

    if dbname == '' or port == '' or dbhost == '' or uname=='':
        raise Exception("Invalid database parameters -- parameters not specified in config file")
            
    ## create connection to db and create necessary tables 
    if verbose:
        print "connecting to database: %s"%dbname
    engine = create_engine('postgresql://%s:%s@%s:%s/%s'%(uname,upass,dbhost,port,dbname),echo=verbose)
    connection = engine.connect()
    Session = sessionmaker(bind=engine)
    session = Session()
    if verbose:
        print 'connected.'

    return session,engine

def read_gene_info_file(lineCount=False,short=False):
    """
    read the essential info from NCBI's gene info file
    """

    config = Configure()
    taxaList = config.log['taxa']
    geneInfoFile = os.path.join(config.log['data'],"gene_info.db")
    geneInfoFid = open(geneInfoFile,'rU')
    header = geneInfoFid.next()
    geneInfo ={}
    totalLines = 0

    for record in geneInfoFid:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]):
            continue

        taxId = record[0]
        ncbiId = record[1]

        if taxId not in taxaList:
            continue

        if lineCount == True:
            totalLines += 1
            continue

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

    config = Configure()
    taxaList = config.log['taxa']
    geneInfoCount = read_gene_info_file(lineCount=True)
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rU')
    reader = csv.reader(idmappingFid,delimiter="\t")
    records = set([])
    totalRecords = 0

    for record in reader:
        totalRecords += 1
        if record[1] == 'ncbi-taxa-id' and record[2] not in taxaList:
            continue
        records.update([record[0]])

    idmappingFid.close()
    print 'debug. totalRecords %s, records %s, geneInfoCount %s'%(totalRecords, len(list(records)), geneInfoCount)
    return len(list(records)),geneInfoCount

def populate_taxon_table(engine):
    """
    given a list of taxon ids populate the taxon table    
    populates all taxa
    """
    
    config = Configure()
    namesFile = os.path.join(config.log['data'],"names.dmp")
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

    config = Configure()
    taxaList = config.log['taxa']
    timeStart = time.time()
    toAdd = []
    totalRecords = 0
    total = geneInfoCount
    wayPoints = [round(int(w)) for w in np.linspace(0,total,20)]
    geneInfoFile = os.path.join(config.log['data'],"gene_info.db")
    geneInfoFid = open(geneInfoFile,'rU')
    header = geneInfoFid.next()
    taxaIdMap = taxa_mapper(session)

    for record in geneInfoFid:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]):
            continue

        taxId = record[0]

        if taxId not in taxaList:
            continue

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
        
        if len(toAdd) >= 200000:
            toRemove = []
            for ta in toAdd:
                if taxaIdMap.has_key(ta['taxa_id']):
                    ta['taxa_id'] = taxaIdMap[ta['taxa_id']]
                else:
                    toRemove.append(ta)

            for ta in toRemove:
                toAdd.remove(ta)
                
            if len(toAdd) > 0:
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

    if len(toAdd) > 0:
        with engine.begin() as connection:
            connection.execute(Gene.__table__.insert().
                               values(toAdd))

    ## clean up
    geneInfoFid.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique genes were added."%totalRecords
    return timeStr,addedStr

def populate_uniprot_table(lineCount,session,engine):
    """
    populate the uniprot table with entries from idmappings
    """

    config = Configure()
    taxaList = config.log['taxa']
    timeStart = time.time()
    totalLines,totalRecords = 0,0
    idmappingFile = get_idmapping_file()
    idmappingFid = open(idmappingFile,'rb')
    reader = csv.reader(idmappingFid,delimiter="\t")
    ac2kbMap,toAdd = {},{}
    wayPoints = [round(int(w)) for w in np.linspace(0,lineCount,20)]

    print("getting mappers...")
    geneIdMap = gene_mapper(session)
    taxonIdMap = taxa_mapper(session)
    print("mappers loaded... %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

    def queue_entries(toAdd,geneIdMap,taxonIdMap,engine):

        toCommit = []

        for uniprotKbEntry, entry in toAdd.iteritems():
            db_gene_id = None
            db_taxa_id = None
            db_gene_taxa_id = None

            ## convert the gene id to a database key (check old names if we cannot find it)
            if entry['gene-id'] == None:
                pass
            elif geneIdMap.has_key(entry['gene-id']):
                db_gene_id = geneIdMap[entry['gene-id']]
            elif not geneIdMap.has_key(entry['gene-id']):
                _geneIds = [re.sub("\s+","",_ncid) for _ncid in entry['gene-id'].split(";")]
                db_gene_id = None
        
                for _gid in _geneIds:
                    if geneIdMap.has_key(_gid):
                        db_gene_id= _gid

            ## convert the taxa id to a database key
            if entry['ncbi-taxa-id'] and taxonIdMap.has_key(entry['ncbi-taxa-id']):
                db_taxa_id = taxonIdMap[entry['ncbi-taxa-id']]

            ## check that the linked gene taxa is the same as the entry taxa
            #if db_gene_id:
            #    db_gene_taxa_id = session.query(Gene).filter_by(id=db_gene_id).first().taxa_id
            #if db_taxa_id and db_gene_id:
            #    if db_taxa_id != db_gene_taxa_id:
            #        print("WARNING: two taxa present in single uniprot entry? %s %s "%(uniprotKbEntry,\
            #                                                                           entry['gene-id']))

            ## if no taxa was provdied use the one assocated with the linked gene
            #if not db_taxa_id:
            #    dgeneQuery = session.query(Gene).filter_by(id=db_gene_id).first()
            #    if dgeneQuery:
            #        db_taxa_id = dgeneQuery.taxa_id 

            ## ensure we are in appropriate taxa
            if entry['ncbi-taxa-id'] not in taxaList:
                continue

            ## ready the uniprot-ac and refseq rows
            entry['uniprot-ac'] = list(entry['uniprot-ac'])
            if len(entry['uniprot-ac']) == 0:
                entry['uniprot-ac'] = None
            elif  len(entry['uniprot-ac']) == 1:
                entry['uniprot-ac'] = entry['uniprot-ac'][0]
            else:
                entry['uniprot-ac'] = ";".join(entry['uniprot-ac'])

            entry['refseq'] = list(entry['refseq'])
            if len(entry['refseq']) == 0:
                entry['refseq'] = None
            elif  len(entry['refseq']) == 1:
                entry['refseq'] = entry['refseq'][0]
            else:
                entry['refseq'] = ";".join(entry['refseq'])

            ## commit to db
            toCommit.append({'uniprot_ac':entry['uniprot-ac'],'uniprot_entry':uniprotKbEntry,
                             'refseq':entry['refseq'],'taxa_id':db_taxa_id,'gene_id':db_gene_id})
        if len(toCommit) > 0:
            with engine.begin() as connection:
                connection.execute(Uniprot.__table__.insert().
                                   values(toCommit))

    ## parse the idmapping file into the db
    for record in reader:

        if len(record) != 3:
            continue

        uniprotKbAc,uniprotKbEntry,ncbiId,refseq,ncbiTaxaId = None,None,None,None,None
        uniprotKbAc = record[0]
        totalLines += 1
        if totalLines in wayPoints:
            print("\t%s / %s"%(totalLines,lineCount))
        
        if record[1] == 'NCBI_TaxID':
            ncbiTaxaId = record[2]
        elif record[1] == 'GeneID':
            ncbiId = record[2]
        elif record[1] == 'UniProtKB-ID':
            uniprotKbEntry = record[2]
            if not ac2kbMap.has_key(uniprotKbAc):
                ac2kbMap[uniprotKbAc] = uniprotKbEntry
        elif record[1] == 'RefSeq':
            refseq = record[2]
        else:
            continue

        ## skip the XXXX-1 like uniprot ac
        if ac2kbMap.has_key(uniprotKbAc) == False:
            continue

        ## get current key
        uniprotKbEntry = ac2kbMap[uniprotKbAc] 

        ## make new entry if necessary
        if uniprotKbEntry and not toAdd.has_key(uniprotKbEntry):

            ## queue entries in blocks
            totalRecords += 1 
        
            if totalRecords % 100000 == 0:
                queue_entries(toAdd,geneIdMap,taxonIdMap,engine)
                toAdd,ac2kbMap = {},{}


            ac2kbMap[uniprotKbAc] = uniprotKbEntry
            toAdd[uniprotKbEntry] = {'ncbi-taxa-id':None,
                                     'gene-id':None,
                                     'uniprot-ac':set([]),
                                     'refseq':set([])}
        
        ## populate uniprot dictionary
        toAdd[uniprotKbEntry]['uniprot-ac'].update([uniprotKbAc])

        if ncbiTaxaId:
            toAdd[uniprotKbEntry]['ncbi-taxa-id'] = ncbiTaxaId
        elif ncbiId:
            toAdd[uniprotKbEntry]['gene-id'] = ncbiId
        elif refseq:
            toAdd[uniprotKbEntry]['refseq'].update([refseq])

    ## queue any remaining
    if len(toAdd.keys()) > 0:
        queue_entries(toAdd,geneIdMap,taxonIdMap,engine)

    ## clean up
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
    config = Configure()
    taxaList = config.log['taxa']
    toAdd = []
    annotationFile = get_annotation_file()
    annotationFid = open(annotationFile,'rU')
    wayPoints = [round(int(w)) for w in np.linspace(0,totalAnnotations,20)]
    annotationCount = 0

    print("...loading mappers")
    termIdMap = goterm_mapper(session)
    taxaIdMap = taxa_mapper(session)
    uniprotIdMap = uniprot_mapper(session)
    print("...populating rows")

    def queue_entry(goId,evidenceCode,pubmedRefs,uniprotId,geneId,taxon,toAdd,mapper,ignoredAnnotations):

        ## remove invalid term ids
        if not termIdMap.has_key(goId):
            queryTerm = session.query(GoTerm).filter_by(alternate_id=goId).first()
            if queryTerm == None:
                return
            go_db_id = queryTerm.id
        else:
            go_db_id = termIdMap[goId]

        ## remove invalid uniprot ids
        if uniprotId and not mapper.has_key(uniprotId):
            return
        if uniprotId:
            uniprot_db_id = mapper[uniprotId]
        else:
            uniprot_db_id = None

        ## remove invalid gene ids
        if geneId and not mapper.has_key(geneId):
            return
        if geneId:
            gene_db_id = mapper[geneId]
        else:
            gene_db_id = None

        ## ignore annotations that have an outdated taxon
        if not taxaIdMap.has_key(taxon):
            ignoredAnnotations += 1
            return

        ## get the taxa foreign key
        taxon_db_id = taxaIdMap[taxon]

        toAdd.append({'go_term_id':go_db_id,'evidence_code':evidenceCode,
                      'pubmed_refs':pubmedRefs,'uniprot_id':uniprot_db_id,
                      'gene_id':gene_db_id,'taxa_id':taxon_db_id})

    ## add annotations from uniprot annotation file
    ignoredAnnotationsUniprot = 0
    print("...getting annotations from gene_association (uniprot)")
    for record in annotationFid:
        record = record[:-1].split("\t")

        ## check that it is a uniprot entry
        if record[0][0] == "!":
            continue
        if record[0] != 'UniProtKB':
            continue
        
        uniprotId = record[1]
        dbObjectSymbol = record[2]
        goId = record[4]
        pubmedRefs = record[5]
        evidenceCode = record[6]
        aspect = record[8]
        uniprotEntry = record[10]
        goTermName = record[11]
        taxon = re.sub("taxon:","",record[12])
        date = record[13]
        assignedBy = record[14]

        if taxon not in taxaList:
            continue

        ## parse the uniprot Entry
        if re.search("\|",uniprotEntry):
            uniprotEntry = re.split("\|",uniprotEntry)[0]

        ## ignore annotations with multiple species
        if re.search("\|",taxon):
            continue

        ## update progress
        annotationCount += 1
        if annotationCount in wayPoints:
            print("\t%s / %s"%(annotationCount,totalAnnotations))

        queue_entry(goId,evidenceCode,pubmedRefs,uniprotEntry,None,taxon,toAdd,
                    uniprotIdMap,ignoredAnnotationsUniprot)

        if len(toAdd) >= 100000: # 100000
            with engine.begin() as connection:
                connection.execute(GoAnnotation.__table__.insert().
                                   values(toAdd))
            toAdd = []

    print('committing final changes...')
    print('ignored annotations after uniprot... %s'%(ignoredAnnotationsUniprot))
    with engine.begin() as connection:
        connection.execute(GoAnnotation.__table__.insert().
                           values(toAdd))

    del uniprotIdMap
    annotationFid.close()
    
    ## add annotations from gene2go
    gene2goFile = get_gene2go_file()
    gene2goFid = open(gene2goFile,'rU')
    ignoredAnnotationsGene = 0 
    print("...getting annotations from gene2go")
    header = gene2goFid.next()
    geneIdMap = gene_mapper(session)
    toAdd = []

    for record in gene2goFid:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]) or len(record) != 8:
            continue
    
        taxon = record[0]
        ncbiId = record[1]
        goId = record[2]
        evidenceCode = record[3]
        qualifier = record[4]
        go_term_description = record[5]
        pubmedRefs = record[6]
        go_aspect = record[7]
        annotationCount += 1

        if taxon not in taxaList:
            continue

        if annotationCount in wayPoints:
            print("\t%s / %s"%(annotationCount,totalAnnotations))

        queue_entry(goId,evidenceCode,pubmedRefs,None,ncbiId,taxon,toAdd,
                    geneIdMap,ignoredAnnotationsGene)

        if len(toAdd) >= 100000: #100000
            with engine.begin() as connection:
                connection.execute(GoAnnotation.__table__.insert().
                                   values(toAdd))
            toAdd = []

    print('ignored annotations after gene2go... %s'%(ignoredAnnotationsGene))
    print('committing final changes...')
    
    with engine.begin() as connection:
        connection.execute(GoAnnotation.__table__.insert().
                           values(toAdd))

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique go annotation entries were added."%annotationCount
    return timeStr,addedStr,(ignoredAnnotationsUniprot,ignoredAnnotationsGene)

def print_db_summary():
    """
    print a summary of rows and tables for the database
    """
    
    print('querying database...')
    config = Configure()
    printstr = ""
    session,engine = db_connect(verbose=False)
    printstr += "\nDATABASE - %s - SUMMARY"%config.log['dbname'] + "\n"
    print("\nDATABASE - %s - SUMMARY"%config.log['dbname'])
    for table in [Taxon,Gene,Uniprot,GoTerm,GoAnnotation]:
        print("There are %s entries in the %s table"%(session.query(table).count(),table.__tablename__))
        printstr += "There are %s entries in the %s table"%(session.query(table).count(),table.__tablename__) + "\n"

    return printstr

def get_idmapping_file():
    """
    check for presence of the annotation file
    raise exception when not found
    return the file path 
    """

    config = Configure()
    dataDir = config.log['data']
    idmappingFile = os.path.join(dataDir,'idmapping.dat.db')
    if os.path.exists(idmappingFile) == False:
        raise Exception("Could not find 'idmapping.dat.db' -- did you run FetchDbData.py?")

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
