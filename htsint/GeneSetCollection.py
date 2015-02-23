#!/usr/bin/python
"""
Container class for gene sets generated from a labels file
"""

__author__ = "Adam Richards"

import os,re,sys,csv
import numpy as np
from htsint.blast import BlastMapper

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

        self.allClusters = np.sort(np.unique(self.labels))
    
    def write(self,blastMap=None, sizeMin=4, sizeMax=100,outFile="genesets.gmt"):
        """
        outFile: specifies the output file path (*.gmt)
        also a *.csv file with gene transcript mapping will be created if a bmap is provided

        blastMap: BlastMap returned after loading summary file in BlastMapper 
        sizeMin: minimum size for a gene set
        sizeMax: maximum size for a gene set
        outFile: outfile path

        """

        ## prepare outfiles
        writer = csv.writer(open(outFile,'w'),delimiter="\t")

        if blastMap:
            outFileMap = re.sub("\.gmt",".csv",outFile)
            writerMap = csv.writer(open(outFileMap,'w'))
            writerMap.writerow(["gene_set","gene_id","mapped_transcripts"])

        ## prepare a blastmap
        if blastMap:
            bm = BlastMapper()
            transcript2gene,gene2transcript = bm.get_dicts(blastMap)

        ## save gene sets to file
        failedThreshold = 0
        realizedGenes = 0
        for _k in self.allClusters:

            numGenes = np.where(self.labels==_k)[0].size
            mappedGenes = set([])
            gsName = "gs-"+str(_k)
            clusterInds = np.where(self.labels==_k)[0]
            clusterGenes = self.genes[clusterInds]
            description = self.get_description(clusterGenes)

            ## map the genes
            if blastMap:
                for gene in clusterGenes:
                    geneTranscripts = set([])
                    if not gene2transcript.has_key(gene):
                        continue
            
                    for taxa,matchedTranscripts in gene2transcript[gene].iteritems():
                        mappedGenes.update(matchedTranscripts)
                        geneTranscripts.update(matchedTranscripts)
                    if blastMap:
                        writerMap.writerow([gsName,gene,";".join(list(geneTranscripts))])


            else:
                mappedGenes = clusterGenes

            mappedGenes = list(mappedGenes)
            realizedGenes += len(mappedGenes)

            if len(mappedGenes) >= sizeMin and len(mappedGenes) <= sizeMax: 
                writer.writerow([gsName,description] + mappedGenes)
            else:
                failedThreshold+=clusterGenes.size

        print("-----------------")
        print("sigma: %s"%self.sigma)
        print("k: %s"%self.k)
        print('Total clusters: %s '%self.allClusters.size)
        percentAccepted = float(self.genes.size-failedThreshold) / float(self.genes.size)
        print("Genes in clusters %s/%s (%s)"%(self.genes.size-failedThreshold,self.genes.size,round(percentAccepted,2)) + "%)")
    
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
        self.k = header1[0].split("=")[1]
        self.sigma = float(header1[1].split("=")[1])
        
        genes,labels = [],[]
        for linja in reader:
            if len(linja) < 2:
                continue
            genes.append(linja[0])
            labels.append(linja[1])

        return np.array(genes),np.array(labels)
        
if __name__ == "__main__":
    print "Running..."
