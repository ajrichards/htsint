#!/usr/bin/env python
"""
ParseBlast.py 
    takes as input an xml results file
    creates a csv summary results file

    can pass file handles for results files 
    instead of file paths if running over 
    multiple files
"""

### make imports
import sys,os,re,time,csv,getopt
from Bio.Blast import NCBIXML

class ParseBlast(object):
    """
    parse the xml files for BLAST results
    """

    def __init__(self,filePath,outDir=".",fhResults=None,fhLog=None):
        """
        outDir    - location of where to put results
        fhResults - results file handle (csv.writer)
        fgLog     - log file handle     (csv.writer)

        """

        ## error checking
        if filePath == None or not re.search("\.xml",filePath):
            raise Exception("invalid results file -  must be of type xml \n%s"%filePath)

        self.filePath = filePath
        self.fid1,self.fid2 = None,None

        ## prepare a log file
        if fhLog == None:
            resultsHomeDir,fileName = os.path.split(filePath)
            self.fid1 = open(os.path.join(outDir,'%s_1.log'%re.sub("\.xml","",fileName)),'w')
            self.logWriter = csv.writer(self.fid1)
        else:
            self.logWriter = fhLog

        ## prepare a summary file
        if fhResults == None:
            summaryFilePath = os.path.join(outDir,'%s_1.csv'%re.sub("\.xml","",fileName))
            self.fid2 = open(summaryFilePath,'w')
            self.resultsWriter = csv.writer(self.fid2)
            self.resultsWriter.writerow(["query","query_length","accession","e-score","bit-score"])
        else:
            self.resultsWriter = fhResults

    def push_out(self,line):
        """
        push a string to both STDOUT and the logfile
        """

        self.logWriter.writerow([line])
        print line

    def run(self):
        """
        run the parser
        currently only saves the best result
        """
        self.push_out(sys.argv[0])
        self.push_out(time.asctime())
        self.push_out("Parsing results file... %s"%self.filePath)

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
                #print query
                #sys.exit()
                #queryContig = re.sub("\s+","",query[0])
                #queryIsogroup = re.sub("\s+","",re.split("\=",query[1])[1])
                #queryLength   = re.sub("\s+","",re.split("\=",query[2])[1])
                #print dir(record.alignments[0])
                queryLength = "nan"
                bestAccession = record.alignments[0].accession
                bestEscore = record.alignments[0].hsps[0].expect
                bestBitScore =  record.alignments[0].hsps[0].score
                print 'blah1', record.alignments[0].title
                print 'blah2', record.alignments[0].hit_id
                print 'blah3', record.alignments[0].hit_def

                self.resultsWriter.writerow([query,bestAccession,bestEscore,bestBitScore])
                hasResults += 1

        self.push_out("total blasted sequences: %s"%totalResults)
        self.push_out("sequences with results : %s"%hasResults)

        ## clean up
        if self.fid1 != None:
            self.fid1.close()
        if self.fid2 != None:
            self.fid2.close()
        print 'complete.'
