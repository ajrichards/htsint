#!/usr/bin/env python
"""
Calcuates the pairwise distances between genes 
"""

__author__ = "Adam Richards"

import os,sys,csv,shutil,cPickle,getopt
import numpy as np
import networkx as nx
from basedir import __basedir__
from multiprocessing import Pool, cpu_count

def mp_worker((source,sink,graphPath)):
    """
    find shortest path length
    """
    G = nx.read_gpickle(graphPath)
    minDistance = 1e8
    if G.has_node(source) and G.has_node(sink):
        (dijkDist, dijkPath) = nx.bidirectional_dijkstra(G,source,sink)
    else:
        dijkDist = None
    
    if dijkDist:
        return [source,sink,dijkDist]

class TermDistances(object):
    """
    a class to represent distances among all gene-ontology terms
    """

    def __init__(self,termsPath,termGraphPath,resultsDir=os.path.join(".","htsint-tmp"),cluster=True):
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
        self.queue = []
        self.baseDir =  os.path.realpath(os.path.dirname(__file__))
            
        ## results dir will be emptied on each run
        self.resultsDir = os.path.realpath(resultsDir)

        ## load the term graph and the terms
        self.G = nx.read_gpickle(self.termGraphPath)
        tmp = open(self.termsPath,'r')
        self.gene2go,self.go2gene = cPickle.load(tmp)
        tmp.close()

        ## variables
        self.terms = self.go2gene.keys()
        self.totalTerms = float(len(self.terms))
        self.totalDistances = ((self.totalTerms*(self.totalTerms+1))/2) - self.totalTerms

    def create_scripts(self,email,name='dist',cpus=1):
        """
        creates bash scripts to be submitted to the cluster queue
        """

        ## empty the results directory
        if os.path.isdir(self.resultsDir):
            shutil.rmtree(self.resultsDir)
        os.mkdir(self.resultsDir)
        
        ## variables
        chunks = cpus
        script = os.path.join(__basedir__,"TermDistances.py")
        stopPoints = np.arange(0,self.totalDistances,
                               int(round(float(self.totalDistances)/float(chunks))))
        if stopPoints[-1] < self.totalDistances:
            stopPoints = np.hstack([stopPoints[1:],np.array([self.totalDistances])])

        print '... submitting %s jobs'%(len(stopPoints))

        ## create scripts
        begin = 0
        for i,chunk in enumerate(range(stopPoints.size)):
            stop = stopPoints[chunk] 

            submitFile = os.path.join(self.resultsDir,"%s-%s.sh"%(name,i))
            submitLog =  os.path.join(self.resultsDir,"%s-%s.log"%(name,i))
            args = " -f %s -l %s -t '%s' -g '%s' -r %s"%(int(begin),int(stop),self.termsPath,
                                                         self.termGraphPath,self.resultsDir)
            f = open(submitFile, 'w')
            f.write("#!/bin/bash\n" +
                    "#$ -S /bin/bash\n" +
                    "#$ -j yes\n" +
                    "#S -M %s\n"%email +
                    "#$ -o %s\n"%submitLog +
                    "/usr/bin/python "+ script + args)

            f.close()
            self.queue.append(submitFile)
            begin = stop

    def submit(self):
        """
        submit the file to the cluster        
        """
    
        if len(self.queue) == 0:
            raise Exception("There are no script to be submitted to the queue\n did you run create_scripts?")

        for queueFileName in self.queue:
            os.system("qsub " + queueFileName)

    def run(self,first=None,last=None):
        """
        Search for all pairwise shortest paths
        """

        ## create a results file 
        outFile = os.path.join(self.resultsDir,"out-%s-%s.npy"%(first,last))
        mat = np.zeros((last-first,3),).astype(str)
        
        if first == None or last == None:
            first = 0
            last = int(self.totalDistances)

        toRun = range(first,last)

        count = -1
        lineCount = -1
        for i,termI in enumerate(self.terms):
            
            if i % 20 == 0:
                print "%s/%s"%(i,len(self.terms))

            for j,termJ in enumerate(self.terms):
                if j >= i:
                    continue
                count += 1
                if count not in toRun:
                    continue
                
                ## find the distance
                distance = self.get_distance(termI,termJ)

                if distance != None:
                    lineCount+=1
                    mat[lineCount,:] = [termI,termJ,distance]

        np.save(outFile,mat)

    def get_distance(self,source,sink):
        """
        find shortest path length
        """

        minDistance = 1e8
        if self.G.has_node(source) and self.G.has_node(sink):
            (dijkDist, dijkPath) = nx.bidirectional_dijkstra(self.G,source,sink)
            return dijkDist
        return None

    def run_with_multiprocessing(self,resultsFilePath,chunkSize=1000,cpus=8):
        """
        method to calculate distances on a single multicore machine
        """

        ## assemble the pairwise terms
        print("loading data")
        pairwiseTerms = []
        for i,termI in enumerate(self.terms):
            for j,termJ in enumerate(self.terms):
                if j >= i:
                    continue
                pairwiseTerms.append([termI,termJ,self.termGraphPath])
                
        print("loaded.")

        start = 0
        stop = chunkSize
        lineCount = 0 
        mat = np.zeros((len(pairwiseTerms),3),).astype(str)
        numChunks = int(np.ceil(float(len(pairwiseTerms) / float(chunkSize))))
        for chunk in range(numChunks):
            p = Pool(cpus)
            _results = p.map_async(mp_worker,pairwiseTerms[start:start+chunkSize])
            results = _results.get()
            p.close()
            p.join()
 
            for row in results:
                if row:
                    lineCount+=1
                    mat[lineCount,:] = row
        
            print("%s"%((start+chunkSize)/float(len(pairwiseTerms)) * 100.0)+"% complete")
            start += chunkSize

        np.save(resultsFilePath,mat)

if __name__ == "__main__":
    ## read in input file      
    if len(sys.argv) < 3:
        print sys.argv[0] + "-f first -l last -q query_file -d database -e evalue -o outdir"
        sys.exit()

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'f:l:t:g:r:')
    except getopt.GetoptError:
        raise Exception(sys.argv[0] + "-f first -l last -t termsPath -g termGraphPath -r resultsDir")
        
    first,last,termsPath,termGraphPath,resultsDir = None,None,None,None,None
    for o,a in optlist:
        if o == '-f':
            first = int(a)
        if o == '-l':
            last  = int(a)
        if o == '-t':
            termsPath = a
        if o == '-g':
            termGraphPath  = a
        if o == '-r':
            resultsDir  = a

    td = TermDistances(termsPath,termGraphPath,resultsDir)
    td.run(first=first,last=last)
