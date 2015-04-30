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
    
    def write(self,blastMap=None, transcriptMin=9, transcriptMax=1000,outFile="genesets.gmt"):
        """
        outFile: specifies the output file path (*.gmt)
        also a *.csv file with gene transcript mapping will be created if a bmap is provided

        blastMap: BlastMap returned after loading summary file in BlastMapper 
        transcriptMin: minimum size for a gene set
        transcriptMax: maximum size for a gene set
        outFile: outfile path

        """

        print("---------------------")
        print('There are %s genes with at least one annotation'%(len(self.gene2go.keys())))
        print('There are %s genes in the labels file'%(len(self.genes)))

        if blastMap:
            bm = BlastMapper()
            bmGenes = bm.print_summary(blastMap)
            gene2transcript = bm.get_gene_dict(blastMap)
            usableGenes = list(set(bmGenes).intersection(set(self.gene2go.keys())))

        if blastMap:
            print('There are %s genes with at least one BLAST hit'%(len(bmGenes)))
            print('There are %s genes that have both a BLAST hit and an annotation'%(len(usableGenes)))
            #print('There are %s genes in clusters with at least one BLAST hits'%(len(set(self.genes).intersection.(set(bmGenes.keys())))))

        ## prepare outfiles
        writer = csv.writer(open(outFile,'w'),delimiter="\t")

        if blastMap:
            outFileMap = re.sub("\.gmt",".csv",outFile)
            writerMap = csv.writer(open(outFileMap,'w'))
            writerMap.writerow(["gene_set","gene_id","mapped_transcripts"])

        ## save gene sets to file
        failedThreshold = 0
        
        for _k in self.allClusters:
            clusterInds = np.where(self.labels==_k)[0]
            clusterGenes = self.genes[clusterInds]
            gsName = "gs-"+str(_k)
            description = self.get_description(clusterGenes)

            ## map the genes
            if blastMap:
                mapped = set([])
                for gene in clusterGenes:
                    if not gene2transcript.has_key(gene):
                        continue
                    geneTranscripts = gene2transcript[gene]
                    geneTranscripts = list(set([re.sub("\.[0-9]$","",g) for g in geneTranscripts]))

                    if blastMap:
                        writerMap.writerow([gsName,gene,";".join(list(geneTranscripts))])
                    mapped.update(geneTranscripts)
                mapped = list(mapped)
            else:
                mapped = clusterGenes

            ### remove non-unique and versioned genes
            #if len(mapped) > 0:
                

            if len(mapped) >= transcriptMin and len(mapped) <= transcriptMax: 
                writer.writerow([gsName,description] + mapped)
            else:
                failedThreshold+=clusterGenes.size

        print("-----------------")
        print("sigma: %s"%self.sigma)
        print("k: %s"%self.k)
        print('Total clusters: %s '%self.allClusters.size)
        percentAccepted = float(self.genes.size-failedThreshold) / float(self.genes.size)
        print("genes pass threshold %s/%s (%s)"%(self.genes.size-failedThreshold,self.genes.size,round(percentAccepted,2)) + "%)")
        
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
