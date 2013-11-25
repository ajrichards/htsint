#!/usr/bin/env python
"""
See CreateDatabase.py before running PrimeDatabase.py

(1) Read the results file produced by ParceNcbiXml
(2) Populate the taxa in the database that have not already
    been added previously.

A logfile of the output is automatically created
"""

### make imports
import sys,os,re,time,csv,getopt
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation
from DatabaseTools import db_connect
from DatabaseTools import populate_taxon_table,populate_gene_table,populate_accession_table,populate_go_tables

class PrimeDatabase(object):
    """
    Look at BLAST results and add any species to database not already present
    """

    def __init__(self,filePath,resultsDir="."):

        ## error checking
        if filePath == None or not re.search("\.csv",filePath):
            raise Exception("invalid results file -  must be of type csv \n%s"%filePath)

        self.filePath = filePath

        ## prepare a log file
        resultsHomeDir,fileName = os.path.split(filePath)
        fid1 = open(os.path.join(resultsDir,'%s_2.log'%re.sub("\_1\.csv","",fileName)),'w')
        self.writer = csv.writer(fid1)

        ## prepare a summary file 
        summaryFilePath = os.path.join(resultsDir,re.sub("\_1","_2",fileName))
        self.fid2 = open(summaryFilePath,'w')
        self.resultsWriter = csv.writer(self.fid2)
        self.resultsWriter.writerow(["query_contig","query_isogroup","query_length","accession","e-score","bit-score","gene_id"])

    def push_out(self,line):
        """
        push a string to both STDOUT and the logfile
        """

        self.writer.writerow([line])
        print line

    def run(self):
        self.push_out(sys.argv[0])
        self.push_out(time.asctime())
        self.push_out("Connecting to the database...")

        ## conect to the database
        session,engine = db_connect(verbose=False)

        ## read the results file to determine the unique gene ids
        self.push_out("Reading the parsed results file...")

        fid3 = open(self.filePath,'r')
        reader = csv.reader(fid3)
        header = reader.next()
        total = 0
        uniqueGeneIds = set([])

        for linja in reader:
            queryContig = linja[0]
            queryIsogroup = linja[1]
            queryLength = linja[2]
            geneID = re.sub("\s","",linja[3])
            eScore = linja[4]
            bitScore = linja[5]
            uniqueGeneIds.update([geneID])
            total += 1
        fid3.close()
        self.push_out("Total gene ids: %s"%total)
        self.push_out("Unique gene ids: %s"%len(uniqueGeneIds))

        geneInfoFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"gene_info.db")
        if os.path.exists(geneInfoFile) == False:
            raise Exception("ERROR: cannot find gene info file")

        geneInfoFID = open(geneInfoFile,'rU')
        header = geneInfoFID.next()
        uniqueTaxa = set([])
        genesFound = set([])

        for record in geneInfoFID:
            record = record.rstrip("\n")
            record = record.split("\t")
    
            if re.search("^\#",record[0]) or len(record) != 15:
                continue

            taxID = record[0]
            ncbi_id = record[1]
            synonyms = record[4]

            if ncbi_id not in uniqueGeneIds:
                continue

            uniqueTaxa.update([taxID])
            genesFound.update([ncbi_id])

        genesNotFound = uniqueGeneIds.difference(genesFound)    
        self.push_out("...Genes found      : %s"%len(genesFound))
        self.push_out("...Unique taxa found: %s"%len(uniqueTaxa))
        self.push_out("...Genes not found  : %s"%len(genesNotFound))

        ## try to match genes not found to old names
        geneHistoryFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"gene_history.db")
        if os.path.exists(geneHistoryFile) == False:
            raise Exception("cannot find gene history file")
            
        geneHistoryFID = open(geneHistoryFile,'rU')
        header = geneHistoryFID.next()
        oldNameMatches = {}
        for record in geneHistoryFID:
            record = record.rstrip("\n")
            record = record.split("\t")
    
            if re.search("^\#",record[0]):
                continue

            tax_id = record[0]
            geneID = record[1]
            discontinuedGeneID = record[2]

            if discontinuedGeneID == '-':
                continue

            if discontinuedGeneID in genesNotFound:
                if geneID == '-':
                    geneID = "no_longer_exists"

                oldNameMatches[discontinuedGeneID] = geneID
                uniqueTaxa.update([tax_id])

        self.push_out("...of the genes not found %s were updated to current names"%len(oldNameMatches.keys()))
        totalFound = len(genesFound) + len(oldNameMatches.values())
        self.push_out("...of the original %s unique genes %s were matched to well-documented ids"%(len(uniqueGeneIds),totalFound))
        self.push_out("The updated number of unique taxa to add to the db is %s"%len(uniqueTaxa))

        ## write the new results summary file with up2date names
        self.push_out("Reading the parsed results file...")

        fid3 = open(self.filePath,'r')
        reader = csv.reader(fid3)
        header = reader.next()

        for linja in reader:
            geneID = re.sub("\s","",linja[3])

            if oldNameMatches.has_key(geneID):
                self.resultsWriter.writerow(linja+[oldNameMatches[geneID]])
            else:
                self.resultsWriter.writerow(linja+[geneID])
        fid3.close()
        self.push_out("New results file has been written.")

        ## clean up 
        self.fid2.close() 

        ## prime database with new taxa
        self.push_out("Getting ready to prime database...")
        taxaList = []
        for taxID in uniqueTaxa:
            query = session.query(Taxon).filter_by(ncbi_id=taxID).first()

            if query == None:
                taxaList.append(taxID)

        self.push_out("Of the %s taxa there are %s new ones"%(len(uniqueTaxa),len(taxaList)))

        if len(taxaList) == 0:
            print "All taxa extracted from results file are already present in database."
            sys.exit()

        ## taxon table
        timeStr,addedStr = populate_taxon_table(taxaList,session)
        self.push_out(timeStr)
        self.push_out(addedStr)

        ## gene table
        timeStr,addedStr = populate_gene_table(taxaList,session)
        self.push_out(timeStr)
        self.push_out(addedStr)

        ## accession table
        timeStr,addedStr = populate_accession_table(taxaList,session)
        self.push_out(timeStr)
        self.push_out(addedStr)

        ## go terms and go annotations tables
        timeStr,addedStr = populate_go_tables(taxaList,session)
        self.push_out(timeStr)
        self.push_out(addedStr)

        self.push_out("DATABASE - SUMMARY")
        self.push_out("There are %s unique taxa "%session.query(Taxon).count())
        self.push_out("There are %s unique genes   "%session.query(Gene).count())
        self.push_out("There are %s unique accessions"%session.query(Accession).count())
        print "\n"

        ## clean up
        self.fid1.close()
        
