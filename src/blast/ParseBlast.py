#!/usr/bin/env python
"""
ParseBlast.py 
    takes as input an xml results file
    creates a csv summary results file

    can pass file handles for results files 
    instead of file paths if running over 
    multiple files


    BLAST against the swissprot database and it returns the protein gi number
    BLAST against the nr database and it returns the protein gi number
    
    The identifier in the results is relevant to the target database.

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
            self.resultsWriter.writerow(["query","identifier","hit-id","e-score","bit-score"])
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
                query = record.query
                identifier = record.alignments[0].accession
                hitID = record.alignments[0].hit_id
                escore = record.alignments[0].hsps[0].expect
                bitScore =  record.alignments[0].hsps[0].score
                
                ## remove any commas and write the results
                self.resultsWriter.writerow([re.sub(",", "", x) for x in [query,str(identifier),hitID]] + [escore,bitScore])
                hasResults += 1

        self.push_out("total blasted sequences: %s"%totalResults)
        self.push_out("sequences with results : %s"%hasResults)

        ## clean up
        if self.fid1 != None:
            self.fid1.close()
        if self.fid2 != None:
            self.fid2.close()
        print 'complete.'
