#!/usr/bin/env python
"""
See CreateDatabase.py before running PrimeDatabase.py

(1) Read the results file produced by ParceNcbiXml
(2) Populate the taxa in the database that have not already
    been added previously.

A logfile of the output is automatically created


This script in addition to priming the database will map 
all of the identifiers (e.g. protein gi) to a gene_id

"""

### make imports
import sys,os,re,time,csv,getopt
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation
from DatabaseTools import db_connect
from DatabaseAppend import DatabaseAppend
from ConversionTools import convert_gene_ids

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
        self.fid1 = open(os.path.join(resultsDir,'%s_2.log'%re.sub("\_1\.csv","",fileName)),'w')
        self.writer = csv.writer(self.fid1)

        ## prepare a summary file 
        summaryFilePath = os.path.join(resultsDir,re.sub("\_1","_2",fileName))
        self.fid2 = open(summaryFilePath,'w')
        self.resultsWriter = csv.writer(self.fid2)
        self.resultsWriter.writerow(["query","hit-identifier","hit-identifier-long","e-score","bit-score","gene_id"])

    def push_out(self,line):
        """
        push a string to both STDOUT and the logfile
        """

        self.writer.writerow([line])
        print line

    def run_protein_gi(self):
        """
        assumes the hit identifier from the BLAST results was a
        protein gi identifier.

        The gene id match is found when possible.
        """

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
        uniqueGIs = set([])
        uniqueGenes = set([])
        giDict = {}

        for linja in reader:
            query = linja[0]
            hitIdentifier = int(linja[1])
            hitIdentifierLong = (linja[2])
            eScore = linja[3]
            bitScore = linja[4]

            ## regular expression to get the gi identifier
            m = re.findall("gi\|[0-9]+",hitIdentifierLong)
            giId = int(m[0][3:])

            uniqueGIs.update([giId])
            total += 1

            giDict[giId] = hitIdentifierLong

        fid3.close()
        self.push_out("Total protein gi identifiers: %s"%total)
        self.push_out("Unique protein gi identifiers: %s"%len(uniqueGIs))

        ## read the gene2accesion file to map to the protein gi to a gene id
        gene2AccFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"gene2accession.db")
        if os.path.exists(gene2AccFile) == False:
            print "ERROR: populate_accession_table() exiting... could not find gene2AccFile"
            return

        gene2AccFID = open(gene2AccFile,'rU')
        header = gene2AccFID.next()
        uniqueTaxa = set([])
        gisFound = set([])
        gi2gene = {}

        for record in gene2AccFID:
            record = record.rstrip("\n")
            record = record.split("\t")

            if re.search("^\#",record[0]) or len(record) != 16:
                continue

            taxID = record[0]
            gene_id = int(record[1])
            protein_gi = record[6]

            if protein_gi == '-':
                continue

            protein_gi = int(protein_gi)

            if protein_gi not in uniqueGIs:
                continue

            uniqueTaxa.update([taxID])
            gisFound.update([protein_gi])
            gi2gene[protein_gi] = gene_id

            #print "FOUND: ", protein_gi
            #print giDict[protein_gi]
            #print "We have one!", protein_gi
            #sys.exit()

        taxaList = convert_gene_ids(gi2gene.values(),'taxid')

        ## write the new results summary file with up2date names
        self.push_out("Reading the parsed results file...")

        fid3 = open(self.filePath,'r')
        reader = csv.reader(fid3)
        header = reader.next()

        for linja in reader:
            hitIdentifierLong = (linja[2])
            m = re.findall("gi\|[0-9]+",hitIdentifierLong)
            gi = int(m[0][3:])
            if gi2gene.has_key(gi) == True:
                self.resultsWriter.writerow(linja+[gi2gene[gi]])
            else:
                self.resultsWriter.writerow(linja+['-'])

        fid3.close()
        self.push_out("New results file has been written.")

        ## use DatabaseAppend to append to database
        self.push_out("Getting ready to prime database...")
        uniqueTaxa = list(set(taxaList))
        self.push_out(list(uniqueTaxa))

        da = DatabaseAppend(list(uniqueTaxa),force=True)
        da.run()

        self.push_out("DATABASE - SUMMARY")
        self.push_out("There are %s unique taxa "%session.query(Taxon).count())
        self.push_out("There are %s unique genes   "%session.query(Gene).count())
        self.push_out("There are %s unique accessions"%session.query(Accession).count())
        print "\n"        

        ## clean up
        self.fid1.close()
        self.fid2.close() 
