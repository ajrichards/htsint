#!/usr/bin/env python
"""
These are helper scripts to populate the database
"""

### make imports
import sys,os,re,time
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation
from config import CONFIG

def check_version():
    """
    ensure that sqlalchemy is at least 0.8
    """
    if float(sqlalchemy.__version__[:-2]) < 0.8:
        print "ERROR: SQLAlchemy version is less than required version"
        sys.exit()

def db_connect(verbose=False):
    """
    generic function to connect to db
    """
    
    check_version()

    ## declare variables
    uname = CONFIG['dbuser']
    upass = CONFIG['dbpass']
    dbhost = CONFIG['dbhost']
    dbname = CONFIG['dbname']
    port = CONFIG['dbport']

    ## create connection to db and create necessary tables 
    print "connecting to database: %s"%dbname
    engine = create_engine('postgres://%s:%s@%s:%s/%s'%(uname,upass,dbhost,port,dbname),echo=verbose)
    connection = engine.connect()
    Session = sessionmaker(bind=engine)
    session = Session()
    print 'connected.'

    return session,engine

def populate_taxon_table(taxonList,session):
    """
    given a list of taxon ids populate the taxon table
    """

    print '\n...populating the taxa table for %s taxa'%len(taxonList)
    namesFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"names.dmp")
    if os.path.exists(namesFile) == False:
        print "ERROR: Cannot find names.dmp... exiting"
        sys.exit()

    namesFID = open(namesFile,'rU')
    taxaCount = 0
    timeStart = time.time()

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
        if taxID in taxonList:
            pass
        elif taxaCount >= 100:
            continue
        
        ## determine if record exists and add common names up until 3
        query = session.query(Taxon).filter_by(ncbi_id=taxID).first()

        ## if record does not exist
        if query == None and scientificName != None:
            taxaCount += 1
            someTaxon = Taxon(taxID,name=scientificName)
            session.add(someTaxon)
        ## if record exists overwrite 
        elif query != None and scientificName != None:
            taxaCount += 1
            query.scientific_name = scientificName
        ## if record exists add a common name 
        elif query != None and commonName != None:
            if  query.common_name_1 == '':
                query.common_name_1 = commonName
            elif commonName not in [query.common_name_1] and query.common_name_2 == '':
                query.common_name_2 = commonName
            elif commonName not in [query.common_name_1,query.common_name_2] and query.common_name_3 == '':
                query.common_name_3 = commonName
        else:
            continue

    session.commit()
    namesFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr =  "...%s unique taxa were added."%taxaCount
    return timeStr, addedStr

def populate_gene_table(taxonList,session):
    """
    given a list of taxon ids populate the gene table
    """
    
    print '\n...populating the genes table for %s taxa'%len(taxonList)
    ## check that all of the taxa are in the taxon table
    for taxID in taxonList:
        query = session.query(Taxon).filter_by(ncbi_id=taxID).first()
        if query == None:
            print "ERROR: populate_gene_table() exiting - not all taxa are present"
            print "...", taxID
            return

    ## update the gene table
    print "...reading original gene_info file"
    print "...this may take some time"
    
    geneInfoFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"gene_info.db")
    if os.path.exists(geneInfoFile) == False:
        print "ERROR: populate_gene_table() exiting... could not find geneInfoFile"
        return

    geneInfoFID = open(geneInfoFile,'rU')
    header = geneInfoFID.next()
    totalRecords = 0
    timeStart = time.time()

    for record in geneInfoFID:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]) or len(record) != 15:
            continue
    
        taxID = record[0]

        if taxID not in taxonList:
            continue

        ncbi_id = record[1]
        symbol = record[2]
        synonyms = record[4]
        chromosome = record[6]
        map_location = record[7]
        description = record[8]  
    
        ## determine if record exists and add common names up until 3
        queryTax = session.query(Taxon).filter_by(ncbi_id=taxID).first()
        queryGene = session.query(Gene).filter_by(ncbi_id=ncbi_id).first()
        taxa_id = queryTax.id
        
        ## if record does not exist
        if queryGene == None:
            totalRecords += 1
            someGene = Gene(ncbi_id,description,symbol,chromosome,map_location,synonyms,taxa_id)
            session.add(someGene)
        ### if record exists overwrite 
        elif queryGene != None:
            taxaCount += 1
            queryGene.ncbi_id = ncbi_id
            queryGene.description = description
            queryGene.symbol = symbol
            queryGene.chromosome = chromosome
            queryGene.map_location = map_location
            queryGene.synonyms = synonyms
            queryGene.taxa_id = taxa_id

    session.commit()
    geneInfoFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique genes were added."%totalRecords
    return timeStr,addedStr

def populate_accession_table(taxonList,session):
    """
    given a list of taxon ids populate the accession table
    """
    
    print '\n...populating the accession table'
    ## check that all of the taxa are in the taxon table
    for taxID in taxonList:
        query = session.query(Taxon).filter_by(ncbi_id=taxID).first()
        if query == None:
            print "ERROR: populate_accession_table() exiting - not all taxa are present"
            print "...", taxID
            return

    ### update the gene table
    print "...reading original gene2accession file"
    print "...this may take a minute"
    
    gene2AccFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"gene2accession.db")
    if os.path.exists(gene2AccFile) == False:
        print "ERROR: populate_accession_table() exiting... could not find gene2AccFile"
        return

    gene2AccFID = open(gene2AccFile,'rU')
    header = gene2AccFID.next()
    totalRecords = 0
    timeStart = time.time()

    for record in gene2AccFID:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]) or len(record) != 16:
            continue
      
        taxID = record[0]
        
        if taxID not in taxonList:
            continue

        ncbi_id= record[1]
        status = record[2]
        rna_nucleo_access_version = record[3]
        rna_nucleo_gi = record[4]
        protein_access_version = record[5]
        protein_gi = record[6]
        genomic_nucleo_access_version = record[7]
        genomic_nucleo_gi = record[8]
        genomic_start = record[9]
        genomic_stop = record[10]
        orientation = record[11] 
        assembly = record[12]

        ## query for the gene key and accession
        queryGene = session.query(Gene).filter_by(ncbi_id=ncbi_id).first()
        gene_ncbi_id = queryGene.id

        ## Record should not exist -- new records will always be appended
        totalRecords += 1
        someAccession = Accession(status,rna_nucleo_gi,protein_gi,genomic_nucleo_gi,genomic_start,
                                  genomic_stop,orientation,assembly,gene_ncbi_id)
        session.add(someAccession)

    session.commit()
    gene2AccFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr =  "...%s unique accession were added."%totalRecords
    return timeStr,addedStr


def populate_go_tables(taxonList,session):
    """
    given a list of taxon ids populate the go tables
    """
    
    print '\n...populating the genes table for %s taxa'%len(taxonList)
    
    ## check that all of the taxa are in the taxon table
    for taxID in taxonList:
        query = session.query(Taxon).filter_by(ncbi_id=taxID).first()
        if query == None:
            print "ERROR: populate_gene_table() exiting - not all taxa are present"
            print "...", taxID
            return

    ## update the gene table
    print "...reading original gene2go file"
    print "...this may take some time"

    gene2goFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"gene2go.db")
    if os.path.exists(gene2goFile) == False:
        print "ERROR: populate_gene_table() exiting... could not find geneInfoFile"
        return

    gene2goFID = open(gene2goFile,'rU')
    header = gene2goFID.next()
    totalTerms,totalAnnotations = 0,0
    timeStart = time.time()

    for record in gene2goFID:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]) or len(record) != 8:
            continue
    
        taxID = record[0]

        if taxID not in taxonList:
            continue

        ncbi_id = record[1]
        go_id = record[2]
        evidence_code = record[3]
        qualifier = record[4]
        go_term_description = record[5]
        pubmed_refs = record[6]
        go_aspect = record[7]
    
        ## determine if record exists and match foreign keys
        queryGene = session.query(Gene).filter_by(ncbi_id=ncbi_id).first()
        gene_id = queryGene.id
        queryTerm = session.query(GoTerm).filter_by(go_id=go_id).first()

        ## if record does not exist add the term
        if queryTerm == None:
            totalTerms += 1
            someTerm = GoTerm(go_id,go_aspect,go_term_description)
            session.add(someTerm)
            queryTerm = session.query(GoTerm).filter_by(go_id=go_id).first()
        
        ## add the annotation
        totalAnnotations+=1
        go_term_id = queryTerm.id
        someAnnotation = GoAnnotation(go_term_id,evidence_code,pubmed_refs,gene_id)
        session.add(someAnnotation)

    session.commit()
    gene2goFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique terms and %s unique annotations were added."%(totalTerms,totalAnnotations)
    return timeStr,addedStr
