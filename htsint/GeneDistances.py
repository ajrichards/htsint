#!/usr/bin/env python
"""
Calcuates the pairwise distances between genes 
"""

__author__ = "Adam Richards"

import os,sys,csv,shutil,cPickle,getopt
import numpy as np
import networkx as nx
from basedir import __basedir__

class GeneDistances(object):
    """
    A generic class
    """

    def __init__(self,termsPath,termGraphPath,resultsDir=os.path.join(".","htsint-tmp")):
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
        self.totalGenes = float(len(self.gene2go.keys()))
        self.totalDistances = ((self.totalGenes*(self.totalGenes+1))/2) - self.totalGenes

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
        script = os.path.join(__basedir__,"GeneDistances.py")
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
        outFile = os.path.join(self.resultsDir,"out-%s-%s.csv"%(first,last))
        outFid = open(outFile,'wa')
        writer = csv.writer(outFid)
        writer.writerow(["i","j","distance"])

        allGenes = sorted(self.gene2go)

        if first == None or last == None:
            first = 0
            last = len(allGenes)

        toRun = range(first,last) 

        count = -1
        for i,geneI in enumerate(allGenes):
            termsI = self.gene2go[geneI]
            
            if i % 20 == 0:
                print "%s/%s"%(i,len(allGenes))

            for j,geneJ in enumerate(allGenes):
                termsJ = self.gene2go[geneJ]
                if j >= i:
                    continue
                count += 1
                if count not in toRun:
                    continue
                
                ## find the distance
                distance = self.get_min_term_dist(termsI,termsJ)
                if distance != None:
                    writer.writerow([geneI,geneJ,distance])

        outfid.close()

    def get_min_term_dist(self,sourceTerms,sinkTerms):
        """
        To get the minimum distance between a source and a sink gene
        we use the minimum distance between their go terms
        return the shorest path between two sets of terms
        """

        minDistance = 1e8
        for source in sourceTerms:
            source = str(source)
            if self.G.has_node(source) == False:
                continue

            for sink in sinkTerms:
                sink = str(sink)
                if source == sink or self.G.has_node(sink) == False:
                    continue
                (dijkDist, dijkPath) = nx.bidirectional_dijkstra(self.G,source,sink)

                if dijkDist < minDistance:
                    minDistance = dijkDist

        if minDistance == 1e8:
            minDistance = None

        return minDistance

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

    gd = GeneDistances(termsPath,termGraphPath,resultsDir)
    gd.run(first=first,last=last)
