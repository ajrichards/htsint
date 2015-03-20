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

    def __init__(self,filePath,outFile=None,fhResults=None,fhLog=None,BLASTDB=None):
        """
        outFile   - location of where to put results
        fhResults - results file handle (csv.writer)
        fgLog     - log file handle     (csv.writer)

        """

        ## error checking
        if filePath == None or not re.search("\.xml|\.outfmt5",filePath):
            raise Exception("invalid results file -  must be of type xml \n%s"%filePath)

        self.filePath = filePath
        self.fid1,self.fid2 = None,None

        ## prepare a log file
        if fhLog == None:
            resultsHomeDir,fileName = os.path.split(filePath)
            self.fid1 = open(os.path.join(".",'%s_parsed.log'%re.sub("\.xml|\.outfmt5","",fileName)),'w')
            self.logWriter = csv.writer(self.fid1)
        else:
            self.logWriter = fhLog

        ## prepare a summary file
        if fhResults == None:
            if outFile == None:
                outFile = os.path.join('.','%s_parsed.csv'%re.sub("\.xml|\.outfmt5","",fileName))
            self.fid2 = open(outFile,'w')
            self.resultsWriter = csv.writer(self.fid2)
            self.resultsWriter.writerow(["query","hit-identifier","hit-identifier-long","e-score","bit-score"])
        else:
            self.resultsWriter = fhResults
        self.outFile = outFile

        ## set the environmental variable if specified  
        if BLASTDB != None:
            os.environ['BLASTDB'] = BLASTDB

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

        hasResults = set([])
        noResults = set([])
        totalHits = 0
        print "parsing... %s"%('0'),
        for record in blast:
            totalHits += 1
            print "\rparsing... %s"%(totalHits),
            if record.alignments:
                query = record.query
                for align in record.alignments:
                    identifier = align.accession
                    hitIdLong = align.title
                    escore = align.hsps[0].expect
                    bitScore =  align.hsps[0].score
                
                    ## remove any commas and write the results
                    self.resultsWriter.writerow([re.sub(",", "", x) for x in [query,str(identifier),hitIdLong]] + [escore,bitScore])
                    hasResults.update([query])
            else:
                noResults.update([record.query])

        hasResults = list(hasResults)
        noResults = list(noResults)
        print("\n")
        self.push_out("total hits: %s"%totalHits)
        self.push_out("sequences with at least one match : %s"%len(hasResults))
        self.push_out("sequences without any matches: %s"%len(noResults))

        ## clean up
        if self.fid1 != None:
            self.fid1.close()
        if self.fid2 != None:
            self.fid2.close()
        print 'complete.'

        return self.outFile
