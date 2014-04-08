#!/usr/bin/env python
"""
These are helper scripts to populate the database
"""

### make imports
import sys,os,re,time,csv
import sqlalchemy
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation

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


def populate_taxon_table(taxonList,session):
    """
    given a list of taxon ids populate the taxon table
    """

    print '\n...populating the taxa table for %s taxa...'%(len(taxonList))
    taxonList = list(set([str(tax) for tax in taxonList]))
    namesFile = os.path.join(CONFIG['data'],"names.dmp")
    if os.path.exists(namesFile) == False:
        print "ERROR: Cannot find names.dmp... exiting"
        sys.exit()

    namesFID = open(namesFile,'rU')
    taxaCount = 0
    timeStart = time.time()
    toAdd = []

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
        else:
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

def populate_gene_table(taxonList,annotations,session):
    """
    use the annotation data to populate all of the genes in the 
    original gene_association.goa_uniprot.gz file.
    
    use the gene_info.gz and gene2unigene
    """

    #taxonList = list(set([str(tax) for tax in taxonList]))
    #print '\n...populating the genes table for %s taxa'%len(taxonList)
    ## check that all of the taxa are in the taxon table
    #for taxID in taxonList:
    #    query = session.query(Taxon).filter_by(ncbi_id=taxID).first()
    #    if query == None:
    #        print "ERROR: populate_gene_table() exiting - not all taxa are present"
    #        print "...", taxID
    #        return

    timeStart = time.time()
    gene2unigeneFile = os.path.join(CONFIG['data'],"gene2unigene")
    geneInfoFile = os.path.join(CONFIG['data'],"gene_info.db")
    
    for fileName in [gene2unigeneFile,geneInfoFile]:
        if not os.path.exists(fileName):
            raise Exception("Could not find the data did you run FetchDbData.py? \n%s"%fileName)

    ## read in the gene2unigene mappings
    gene2unigeneFid = open(gene2unigeneFile,'rU')
    unigene2gene = {}
    for record in gene2unigeneFid:
        record = record[:-1].split("\t")
        if record[0][0] == "#":
            continue
        unigene2gene[record[1]] = record[0]
    gene2unigeneFid.close()

    debug = 0 
    print unigene2gene.keys()[:10]
    for uniprotId, items in annotations.iteritems():
        debug += 1

        geneId = None
        for unigeneName in items['names']:
            if unigene2gene.has_key(unigeneName):
                geneId = unigene2gene[unigeneName]

        print uniprotId, geneId, items.keys()



        if debug == 10:
            sys.exit()

        print "\t", items['names']

        #someGene = Gene(ncbi_id,description,symbol,chromosome,map_location,synonyms,taxa_id)
        

    ## read in the gene info file
    #geneInfoDict = {}
    #geneInfoFid = open(geneInfoFile,'rU')
    #header = geneInfoFid.next()
    #totalRecords = 0
    #
    #for record in geneInfoFid:
    #    record = record.rstrip("\n")
    #    record = record.split("\t")
    #
    #    if re.search("^\#",record[0]) or len(record) != 15:
    #        continue
    #
    #    taxID = record[0]
    #    if taxID not in taxonList:
    #        continue
    #    
    #print 'done'
        


    sys.exit()

    session.add_all(toAdd)
    session.commit()
    geneInfoFid.close()



    '''
    ## update the gene table
    print "...reading original gene_info file"
    print "...this may take some time"
    taxonList = list(set(taxonList))
    geneInfoFile = os.path.join(CONFIG['data'],"gene_info.db")
    if os.path.exists(geneInfoFile) == False:
        print "ERROR: populate_gene_table() exiting... could not find geneInfoFile"
        return

    toAdd = []
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
            toAdd.append(someGene)
        ### if record exists overwrite 
        elif queryGene != None:
            queryGene.ncbi_id = ncbi_id
            queryGene.description = description
            queryGene.symbol = symbol
            queryGene.chromosome = chromosome
            queryGene.map_location = map_location
            queryGene.synonyms = synonyms
            queryGene.taxa_id = taxa_id

        if len(toAdd) > 10000:
            session.add_all(toAdd)
            toAdd = []

    session.add_all(toAdd)
    session.commit()
    geneInfoFID.close()
    '''

    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique genes were added."%totalRecords
    return timeStr,addedStr
    

def populate_accession_table(taxonList,session):
    """
    given a list of taxon ids populate the accession table
    """
    
    print '\n...populating the accession table'
    taxonList = list(set([str(tax) for tax in taxonList]))

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
    
    gene2AccFile = os.path.join(CONFIG['data'],"gene2accession.db")
    if os.path.exists(gene2AccFile) == False:
        print "ERROR: populate_accession_table() exiting... could not find gene2AccFile"
        return

    gene2AccFID = open(gene2AccFile,'rU')
    header = gene2AccFID.next()
    totalRecords = 0
    timeStart = time.time()
    toAdd = []

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
        toAdd.append(someAccession)

        if len(toAdd) > 1000:
            session.add_all(toAdd)
            toAdd = []

    session.add_all(toAdd)
    session.commit()
    gene2AccFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr =  "...%s unique accession were added."%totalRecords
    return timeStr,addedStr


def populate_go_tables(taxonList,session):
    """
    given a list of taxon ids populate the go tables
    """

    print '\n...populating the go term and annotation tables for taxa'
    taxonList = list(set([str(tax) for tax in taxonList]))

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
    gene2goFile = os.path.join(CONFIG['data'],"gene2go.db")
    if os.path.exists(gene2goFile) == False:
        print "ERROR: populate_gene_table() exiting... could not find geneInfoFile"
        return

    gene2goFID = open(gene2goFile,'rU')
    header = gene2goFID.next()
    totalTerms,totalAnnotations = 0,0
    timeStart = time.time()
    toAddAnnotations = []

    for record in gene2goFID:
        record = record.rstrip("\n")
        record = record.split("\t")

        if re.search("^\#",record[0]) or len(record) != 8:
            continue
    
        taxID = record[0]

        if taxID not in taxonList:
            continue

        ## determine if record exists and add common names up until 3
        queryTax = session.query(Taxon).filter_by(ncbi_id=taxID).first()
        taxa_id = queryTax.id
        
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
        someAnnotation = GoAnnotation(go_term_id,evidence_code,pubmed_refs,gene_id,taxa_id)
        toAddAnnotations.append(someAnnotation)

        if len(toAddAnnotations) > 1000:
            session.add_all(toAddAnnotations)
            toAddAnnotations = []

    session.add_all(toAddAnnotations)
    session.commit()
    gene2goFID.close()
    timeStr = "...total time taken: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    addedStr = "...%s unique terms and %s unique annotations were added."%(totalTerms,totalAnnotations)
    return timeStr,addedStr







def print_go_summary(outfile=os.path.join(".","go_summary.csv")):
    """
    print a summary about the organisms with gene ontology information
    """

    session,engine = db_connect()
    goTaxa = []#get_all_go_taxa()
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
