#!/usr/bin/python
"""
The output from ParallelBlast.py is a large numbers of files.
This is a convenience class for parsing these results as 
we would normally with ParseBlast.py.
"""
import shutil,os,sys,re,csv
import numpy as np
from Bio import SeqIO
from htsint import __basedir__
from htsint.blast import ParseBlast

class ParseParallelBlast(object):
    
    def __init__(self,queryFile,resultsDir=os.path.join(".","cluster")):
        """
        Constructor
        queryFile - is a fasta file of sequences
        resultsDir - is a directory that is purged and filled
        """

        ## variables
        self.queryFile = os.path.realpath(queryFile)
        self.resultsDir = os.path.realpath(resultsDir)
        self.evalue = 0.05

        ## ensure the results directory is present and contains files
        if os.path.isdir(resultsDir) == False or len(os.listdir(resultsDir)) == 0:
            raise Exception("The results directory does not exist")

    def parse_file(self,resultsFilePath):
        """
        use ParseBlast to parse a given xml formatted file
        """

        print '...parsing', resultsFilePath
        parser = ParseBlast(resultsFilePath,fhLog=self.logWriter,fhResults=self.resultsWriter)
        parser.run()

    def run(self,chunks):
        """
        runs ParseBlast on each file and combines the results
        """

        ## variables
        queryFileName = os.path.split(self.queryFile)[-1]
        resultsNameBase = re.sub("\.\w+","",queryFileName,flags=re.IGNORECASE)

        ## prepare a log outfile
        self.fid1 = open(os.path.join(self.resultsDir,resultsNameBase+'_1.log'),'w')
        self.logWriter = csv.writer(self.fid1)
        
        ## prepare a results outfile
        self.fid2 = open(os.path.join(self.resultsDir,resultsNameBase+'_1.csv'),'w')
        self.resultsWriter = csv.writer(self.fid2)
        self.resultsWriter.writerow(["query","hit-identifier","hit-identifier-long","e-score","bit-score"])
        
        ## get the start and stop points
        handleIn = open(self.queryFile, "rU")
        total = 0
        for record in SeqIO.parse(handleIn,"fasta") :
            total += 1
            
        handleIn.close()
        stopPoints = np.arange(0,total,int(round(float(total)/float(chunks))))
        if stopPoints[-1] < total:
            stopPoints = np.hstack([stopPoints[1:],np.array([total])])
     
        ## get the name of the results files
        begin = 0
        for i,chunk in enumerate(range(stopPoints.size)):
            stop = stopPoints[chunk] 
            resultsFile = resultsNameBase + "-%s-%s.xml"%(begin,stop)
            resultsFilePath = os.path.join(self.resultsDir,resultsFile)
            if os.path.exists(resultsFilePath)  == False:
                raise Exception("Cannot find one of the chunked results files -- did you imput the correct number of chunks?")
            self.parse_file(resultsFilePath)
            begin = stop

        ## cleanup 
        self.fid1.close()
        self.fid2.close()
        print "\nParsing completed."
        print "Log file: ", os.path.join(self.resultsDir,resultsNameBase+'_1.log')
        print "Results : ", os.path.join(self.resultsDir,resultsNameBase+'_1.csv')                                 
