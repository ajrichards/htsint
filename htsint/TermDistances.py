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


def get_distance_mp(args):
    """
    find shortest path length (same as class function except abstracted out to work with multiprocessing)
    """

    source,sink,G = args
    minDistance = 1e8
    if G.has_node(source) and G.has_node(sink):
        (dijkDist, dijkPath) = nx.bidirectional_dijkstra(G,source,sink)
        return dijkDist
    else:
        return None

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
        outFile = os.path.join(self.resultsDir,"out-%s-%s.csv"%(first,last))
        outFid = open(outFile,'wa')
        writer = csv.writer(outFid)
        writer.writerow(["i","j","distance"])

        if first == None or last == None:
            first = 0
            last = int(self.totalDistances)

        toRun = range(first,last)

        count = -1

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
                    writer.writerow([termI,termJ,distance])

        outFid.close()

    def get_distance(self,source,sink):
        """
        find shortest path length
        """

        minDistance = 1e8
        if self.G.has_node(source) and self.G.has_node(sink):
            (dijkDist, dijkPath) = nx.bidirectional_dijkstra(self.G,source,sink)
            return dijkDist
        return None

    def run_with_multiprocessing(self,resultsFilePath):

        ## create a results file )
        outFid = open(resultsFilePath,'wa')
        writer = csv.writer(outFid)
        writer.writerow(["i","j","distance"])

        ## assemble the pairwise terms
        pairwiseTerms = []
        for i,termI in enumerate(self.terms):
            for j,termJ in enumerate(self.terms):
                if j >= i:
                    continue
                pairwiseTerms.append((termI,termJ,self.G))
                
        pairwiseTerms = pairwiseTerms[:100]

        ## run using multiprocessing
        po = Pool(processes=cpu_count()-1)
        _results = po.map_async(get_distance_mp,pairwiseTerms)
        results =  _results.get()
        
        for pt in pairwiseTerms:
            if results[i] != None:
                writer.writerow([pt[0],pt[1],results[i]])

        #po = Pool(processes=cpu_count()-1)
        #_results = po.map_async(self.run(0,self.totalDistancesgreat_circle,(mat[i,:] for i in range(mat.shape[0])))
        #results =  _results.get()

        outFid.close()

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
