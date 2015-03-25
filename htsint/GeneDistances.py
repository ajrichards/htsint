#!/usr/bin/env python
"""
Calcuates the pairwise distances between genes 
"""

__author__ = "Adam Richards"

import os,sys,csv,shutil,cPickle,getopt
import numpy as np
from basedir import __basedir__

class GeneDistances(object):
    """
    A generic class to handle calculate gene distances using term distances
    """

    def __init__(self,termsPath,termDistancesPath,outFile=None):
        """
        Constructor
        """

        ## error checking
        for path in [termsPath]:
            if os.path.exists(path) == False:
                raise Exception("Cannot find specified path\n%s"%path)
        
        ## variables
        self.termsPath = os.path.realpath(termsPath)
        self.queue = []
        self.baseDir =  os.path.realpath(os.path.dirname(__file__))
            
        if outFile == None:
            self.outFile = "gdistances.csv"
        else:
            self.outFile = outFile

        ## load the graph, the terms and the term distances
        tmp = open(self.termsPath,'r')
        self.gene2go,self.go2gene = cPickle.load(tmp)
        tmp.close()

        ## read in the term distances
        termDist = {}
        mat = np.load(termDistancesPath)
        print mat
        print mat.shape


        for i in range(mat.shape[0]):
            
            linja = mat[i,:]

            if not termDist.has_key(linja[0]):
                termDist[linja[0]] = {}
            if not termDist[linja[0]].has_key(linja[1]):
                termDist[linja[0]][linja[1]] = float(linja[2])

        self.termDist = termDist

        ## variables
        self.genes = self.gene2go.keys()
        self.totalGenes = float(len(self.genes))
        self.totalDistances = ((self.totalGenes*(self.totalGenes+1))/2) - self.totalGenes

    def get_min_term_dist(self,sourceTerms,sinkTerms):
        """
        To get the minimum distance between a source and a sink gene
        we use the minimum distance between their go terms
        return the shorest path between two sets of terms
        """

        minDistance = 1e8
        td = 1e8
        for source in sourceTerms:            
            for sink in sinkTerms:
                if source == sink:
                    continue
            
                ## get dist
                if self.termDist.has_key(source) and self.termDist[source].has_key(sink):
                    td = self.termDist[source][sink]
                elif self.termDist.has_key(sink) and self.termDist[sink].has_key(source):
                    td = self.termDist[sink][source]

                if td < minDistance:
                    minDistance = td

        if minDistance == 1e8:
            minDistance = None

        return minDistance


    def run(self,first=None,last=None):
        """
        Search for all pairwise shortest paths
        """

        ## create a results file 
        outFid = open(self.outFile,'wa')
        writer = csv.writer(outFid)
        writer.writerow(["i","j","distance"])

        for i,geneI in enumerate(self.genes):
            termsI = self.gene2go[geneI]
            
            if i % 20 == 0:
                print "%s/%s"%(i,int(self.totalGenes))

            for j,geneJ in enumerate(self.genes):
                termsJ = self.gene2go[geneJ]
                if j >= i:
                    continue
                
                ## find the distance
                distance = self.get_min_term_dist(termsI,termsJ)
                if distance != None:
                    writer.writerow([geneI,geneJ,distance])

        outFid.close()
