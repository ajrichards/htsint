#!/usr/bin/env python
"""
Assembles the files containing the pairwise distances
that are produced with 'GeneDistances.py' 
"""

__author__ = "Adam Richards"

import os,sys,csv,shutil,re,cPickle,getopt
import numpy as np
import networkx as nx
from basedir import __basedir__

class AssembleDistances(object):
    """
    A generic class
    """

    def __init__(self,termsPath,termGraphPath,resultsDir=os.path.join(".","htsint-tmp"),
                 resultsPath=os.path.join(".","assembled-results.csv")):
        """
        Constructor
        
        if cpus > 1 then the script assumes we are in a cluster envrionment
        Note that the results dir will be cleaned before each run
        """

        ## error checking
        for path in [termGraphPath, termsPath]:
            if os.path.exists(path) == False:
                raise Exception("Cannot find specified path\n%s"%path)
        
        ## variables
        self.termsPath = os.path.realpath(termsPath)
        self.termGraphPath = os.path.realpath(termGraphPath)
        self.baseDir =  os.path.realpath(os.path.dirname(__file__))
            
        ## results dir must contain results
        self.resultsDir = os.path.realpath(resultsDir)

        if re.search("\.csv",resultsPath):
            ## create a files to assemple the results in
            self.writer = csv.writer(open(resultsPath,'w'))
            self.writer.writerow(["i","j","distance"])
        else:
            self.resultsPath = resultsPath

        ## load the term graph and the terms
        self.G = nx.read_gpickle(self.termGraphPath)
        tmp = open(self.termsPath,'r')
        self.gene2go,self.go2gene = cPickle.load(tmp)
        tmp.close()

        ## variables
        self.terms = self.go2gene.keys()
        self.totalTerms = float(len(self.terms))
        self.totalDistances = ((self.totalTerms*(self.totalTerms+1))/2) - self.totalTerms
        self.appendedDistances = 0

    def run(self,name='dist',cpus=1):
        """
        creates bash scripts to be submitted to the cluster queue
        """

        ## check that there are results in the results dir
        if not os.path.isdir(self.resultsDir) or len(os.listdir(self.resultsDir)) == 0:
            raise Exception("There are no results in the results dir \n%s"%self.resultsDir)
                           
        ## variables
        chunks = cpus
        stopPoints = np.arange(0,self.totalDistances,
                               int(round(float(self.totalDistances)/float(chunks))))
        if stopPoints[-1] < self.totalDistances:
            stopPoints = np.hstack([stopPoints[1:],np.array([self.totalDistances])])

        print '...assembling results - %s jobs with %s distances'%(len(stopPoints),self.totalDistances)

        ## create scripts
        mat = np.zeros((self.totalDistances,3),).astype(str)
        begin = 0
        for i,chunk in enumerate(range(stopPoints.size)):
            stop = stopPoints[chunk] 
            submitFile = os.path.join(self.resultsDir,"%s-%s.sh"%(name,i))
            submitLog =  os.path.join(self.resultsDir,"%s-%s.log"%(name,i))       
            
            ## determine if we are using csv or npy
            outFile1 = os.path.join(self.resultsDir,"out-%s-%s.csv"%(int(begin),int(stop)))
            outFile2 = os.path.join(self.resultsDir,"out-%s-%s.npy"%(int(begin),int(stop)))

            if os.path.exists(outFile1):
                self.append_file(outFile1)
            elif os.path.exists(outFile2):
                rows = np.load(outFile2)
                print('adding %s rows...'%(int(stop)-int(begin)))
                self.appendedDistances += (int(stop)-int(begin))
                mat[int(begin):int(stop),:] = rows
            else:
                raise Exception("Results file does not exist -- is the number of cpus correct?\n...%s"%outFile)

            begin = stop

        ## print total
        print("saving...")
        np.save(self.resultsPath,mat)
        print "%s/%s distances appended"%(self.appendedDistances,self.totalDistances)

    def append_file(self,outFile):
        """
        appends a file to the results file
        """
        
        print("...appending %s"%outFile)
        fid = open(outFile,'r')
        reader = csv.reader(fid)
        header = reader.next()

        for linja in reader:
            self.writer.writerow(linja)
            self.appendedDistances += 1

        fid.close()
