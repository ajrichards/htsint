#!/usr/bin/python
"""
Container class for gene sets generated from a labels file
"""

__author__ = "Adam Richards"

import os,sys,csv
import numpy as np
from htsint.blast import BlastMapper

#import sys,os,getopt,csv,cPickle,shutil
#import matplotlib.pyplot as plt
#import numpy as np
#from htsint.database import Uniprot,db_connect


## database 
#from htsint.database import db_connect,Gene,Taxon,Uniprot
#session, engine = db_connect()
#conn = engine.connect()


class GeneSetCollection(object):
    """
    gene set class

    """

    def __init__(self, labelsFile, gene2go):
        """
        Constructor

        labelsFile: path to a csv with the following format
 
            k=300,sigma=0.5
            gene,label
            12345,0
            12346,1
            12347,2
            12349,0
            ...
        
        gene2go: a dictionary with the following format

            gene2go = {'100158436': ['GO:0008277', 'GO:0071108', 'GO:0016579',
                       '734137': ['GO:0006508'], 
                       '398019': ['GO:0007281']}

        """

        ## error check
        if not os.path.exists(labelsFile):
            raise Exception("cannot find labels file")

        self.genes,self.labels = self.read_labels(labelsFile)
        print("loaded file with %s genes and %s clusters"%(self.genes.size,np.unique(self.labels).size))
        
        if type(gene2go) != type({}):
            raise Exception("invalid gene2go dictionary")
        self.gene2go = gene2go

        ## get cluster sizes
        clusters = []
        for _k in np.arange(self.k):
            clusters.append(len(np.where(self.labels==_k)[0]))
        self.clusters = np.array(clusters)

    
    def write(self,blastMap=None, sizeMin=4, sizeMax=100,outFile="genesets.gmt"):
        """
        outFile: specifies the output file path (*.gmt)

        blastMap: BlastMap returned after loading summary file in BlastMapper 
        sizeMin: minimum size for a gene set
        sizeMax: maximum size for a gene set
        outFile: outfile path

        """

        ## prepare outfile
        writer = csv.writer(open(outFile,'w'),delimiter="\t")

        ## prepare a blastmap
        if blastMap:
            bm = BlastMapper()
            transcript2gene,gene2transcript = bm.get_dicts(blastMap)

        ## get percentage
        #tooSmall = np.where(self.clusters < sizeMin)[0]
        #tooLarge = np.where(self.clusters > sizeMax)[0]
        #tooSmallGenes = float(self.clusters[tooSmall].sum())
        #tooLargeGenes = float(self.clusters[tooLarge].sum())
        #numer = self.clusters.sum() - tooSmallGenes - tooLargeGenes
        #percentAccepted = ( numer / self.clusters.sum()) * 100.0 

        ## save gene sets to file
        failedThreshold = 0
        realizedGenes = 0
        for _k in np.arange(self.k):

            numGenes = self.clusters[_k]
            mappedGenes = set([])
            gsName = "gs-"+str(int(_k))
            clusterInds = np.where(self.labels==_k)[0]
            clusterGenes = self.genes[clusterInds]
            description = self.get_description(clusterGenes)

            ## map the genes
            if blastMap:
                for gene in clusterGenes:
                    if not gene2transcript.has_key(gene):
                        continue
            
                    for taxa,matchedTranscripts in gene2transcript[gene].iteritems():
                        mappedGenes.update(matchedTranscripts)
            else:
                mappedGenes = clusterGenes

            mappedGenes = list(mappedGenes)
            realizedGenes += len(mappedGenes)

            if numGenes >= sizeMin and numGenes <= sizeMax: 
                writer.writerow([gsName,description] + mappedGenes)
            else:
                failedThreshold+=len(mappedGenes)

        print("-----------------")
        print("sigma: %s"%self.sigma)
        print("k: %s"%self.k)
        print('Total clusters: %s '%self.clusters.size)
        percentAccepted = float(realizedGenes-failedThreshold) / float(realizedGenes)
        print("Genes in clusters %s/%s (%s"%(realizedGenes-failedThreshold,realizedGenes,round(percentAccepted,2)) + "%)")
    
    def get_description(self,geneList):
        """
        function the return description
        """

        terms = []
        for gene in geneList:
            terms.extend(self.gene2go[gene])

        terms = np.array(terms)
        uniqueTerms = np.sort(np.unique(terms))
        counts = []
        for t in uniqueTerms:
            counts.append(np.where(terms==t)[0].size)
        count = np.array(counts)
        rankedCounts = np.argsort(counts)[::-1]
        topThree = uniqueTerms[rankedCounts[:3]].tolist()

        return ';'.join(topThree)

    def read_labels(self, labelsFile):
        """
        read the labels file
        """

        ## read in the genes and labels
        fid = open(labelsFile,'r')
        reader = csv.reader(fid)
        header1 = reader.next()
        header2 = reader.next()
        self.k = int(header1[0].split("=")[1])
        self.sigma = float(header1[1].split("=")[1])
        
        genes,labels = [],[]
        for linja in reader:
            genes.append(linja[0])
            labels.append(int(linja[1]))

        return np.array(genes),np.array(labels)
        
if __name__ == "__main__":
    print "Running..."
