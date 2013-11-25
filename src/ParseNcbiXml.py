#!/usr/bin/env python
"""
ParseNcbiXml.py 
    takes as input an xml results file
    creates a csv summary results file
"""

### make imports
import sys,os,re,time,csv,getopt
from Bio.Blast import NCBIXML


class ParseNcbiXml(object):
    """
    parse the xml files for BLAST results
    """

    def __init__(self,filePath,resultsDir="."):

        ## error checking
        if filePath == None or not re.search("\.xml",filePath):
            raise Exception("invalid results file -  must be of type xml \n%s"%filePath)

        self.filePath = filePath

        ## prepare a log file
        resultsHomeDir,fileName = os.path.split(filePath)
        self.fid1 = open(os.path.join(resultsDir,'%s_1.log'%re.sub("\.xml","",fileName)),'w')
        self.writer = csv.writer(self.fid1)

        ## prepare a summary file
        summaryFilePath = os.path.join(resultsDir,'%s_1.csv'%re.sub("\.xml","",fileName))
        self.fid2 = open(summaryFilePath,'w')
        self.resultsWriter = csv.writer(self.fid2)
        self.resultsWriter.writerow(["query_contig","query_isogroup","query_length","accession","e-score","bit-score"])

    def push_out(self,line):
        """
        push a string to both STDOUT and the logfile
        """

        self.writer.writerow([line])
        print line

    def run(self):
        """
        run the parser
        """
        self.push_out(sys.argv[0])
        self.push_out(time.asctime())
        self.push_out("Parsing results file...")

        result_handle = open(self.filePath)
        blast = NCBIXML.parse(result_handle)

        hasResults = 0
        totalResults = 0
        print "parsing... %s"%('0'),
        for record in blast:
            totalResults += 1
            print "\rparsing... %s"%(totalResults),
        if record.alignments:
            query =  re.split("\s+",record.query)
        
            queryContig = re.sub("\s+","",query[0])
            queryIsogroup = re.sub("\s+","",re.split("\=",query[1])[1])
            queryLength   = re.sub("\s+","",re.split("\=",query[2])[1])
            bestAccession = record.alignments[0].accession
            bestEscore = record.alignments[0].hsps[0].expect
            bestBitScore =  record.alignments[0].hsps[0].score
            self.resultsWriter.writerow([queryContig,queryIsogroup,queryLength,bestAccession,bestEscore,bestBitScore])
            hasResults += 1

        ## debug 
        #if hasResults >= 10:
        #    sys.exit()

        print "\rparsing... %s"%(totalResults)
        self.push_out("total blasted sequences: %s"%totalResults)
        self.push_out("sequences with results : %s"%hasResults)

        ## clean up
        self.fid1.close()
        self.fid2.close()
        print 'complete.'
