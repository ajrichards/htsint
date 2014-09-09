#!/usr/bin/python
"""
break a blast problem into 'chunks' and send the chunks to the cluster
using the pairwise distance problem here

'chunks' may be more than the number of available cpus

currently this is only working with blastx

"""

import shutil,os,sys
import numpy as np
from Bio import SeqIO
from htsint import __basedir__
from htsint.blast import Blast

class ParallelBlast(object):
    
    def __init__(self,queryFile,database,BLASTDB=None,resultsDir=os.path.join(".","cluster"),cmd='blastx'):
        """ 
        Constructor
        queryFile - is a fasta file of sequences
        resultsDir - is a directory that is purged and filled
        """

        ## variables
        self.queryFile = os.path.realpath(queryFile)
        self.database = database
        self.resultsDir = os.path.realpath(resultsDir)
        self.submitFileNames = []
        self.evalue = 0.05
        self.cmd = cmd

        ## ensure the results directory present and clean
        for dirName in [self.resultsDir]:
            if os.path.isdir(dirName):
                shutil.rmtree(dirName)
            os.mkdir(dirName)

    def create_scripts(self,chunks,email,name='blst'):
        """
        chunks - the number of times pieces to break the queue file into
        email  - tells the queue the user email
        name   - a short name (to be seen in the queue)

        """

        ## variables
        script = os.path.join(__basedir__,"blast","Blast.py")

        handleIn = open(self.queryFile, "rU")
        total = 0
        for record in SeqIO.parse(handleIn,"fasta") :
            total += 1
            
        handleIn.close()

        ## break the graph into chunks
        stopPoints = np.arange(0,total,int(round(float(total)/float(chunks))))
        if stopPoints[-1] < total:
            stopPoints = np.hstack([stopPoints[1:],np.array([total])])
        
        print '... submitting %s jobs'%(len(stopPoints))
             
        begin = 0
        for i,chunk in enumerate(range(stopPoints.size)):
            stop = stopPoints[chunk] 
    
            ## make a file for each chunk
            submitFile = os.path.join(self.resultsDir,"%s-%s.sh"%(name,i))
            submitLog =  os.path.join(self.resultsDir,"%s-%s.log"%(name,i))
            args = " -f %s -l %s -q %s -d %s -e %s -o %s -c %s"%(begin,stop,self.queryFile,self.database,
                                                                 self.evalue,self.resultsDir,self.cmd)
            f = open(submitFile, 'w')
            f.write("#!/bin/bash\n" + 
                    "#$ -S /bin/bash\n" +
                    "#$ -j yes\n" +
                    "#S -M %s\n"%email +
                    "#$ -o %s\n"%submitLog + 
                    "/usr/bin/python "+ script + args)
            f.close()
            begin = stop
            self.submitFileNames.append(submitFile)

    def submit(self):
        """
        submit the scripts via qsub
        """
        
        if len(self.submitFileNames) == 0:
            raise Exception("There are no script to be submitted to the queue\n did you run create_scripts?")

        for submitFile in self.submitFileNames:
            os.system("qsub " + submitFile)
