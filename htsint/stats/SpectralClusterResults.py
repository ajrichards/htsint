#!/usr/bin/env python
"""
A generic template
"""

import os,csv
import numpy as np

__author__ = "Adam Richards"

class SpectralClusterResults(object):
    """
    A class to handle spectral clustering results
    """

    def __init__(self,silvalsFile,clustersFile):
        """
        Constructor
        """

        ## error checking
        for filePath in[silvalsFile,clustersFile]:
            if not os.path.exists(filePath):
                raise Exception("could not find file: %s"%filePath)

        self.clusters = self.load_clusters_file(clustersFile)
        self.kRange,self.sigRange,self.silvals = self.load_silvals_file(silvalsFile)

    def load_clusters_file(self,clusterFile):
        """
        load the clusters file
        k,sigma,clustid1,clustid2,....clustid_max

        """
        
        fid = open(clusterFile,'r')
        reader = csv.reader(fid)
        header = reader.next()
        clusterIds = np.array(header[2:])
        results = {}

        for linja in reader:
            k = str(int(linja[0]))
            sigma = str(round(float(linja[1]),4))
            if not results.has_key(k):
                results[k] = {}

            results[k][sigma] = np.array([float(i) for i in linja[2:]])

        fid.close()        
        return results

    def load_silvals_file(self,silvalsFile):
        """
        load the clusters file
        k,sigma,clustid1,clustid2,....clustid_max

        """
        
        fid = open(silvalsFile,'r')
        reader = csv.reader(fid)
        header = reader.next()
        kRange = set([])
        sigRange = set([])
        for linja in reader:
            k = int(linja[0])
            sigma = round(float(linja[1]),4)
            kRange.update([k])
            sigRange.update([sigma])
        fid.close()

        kRange = np.sort(np.array(list(kRange)))
        sigRange = np.sort(np.array(list(sigRange)))

        ## create matrix with k as rows and sigma as columns
        resultsMat = np.zeros((kRange.size,sigRange.size),)

        fid = open(silvalsFile,'r')
        reader = csv.reader(fid)
        header = reader.next()

        for linja in reader:
            k = int(linja[0])
            sigma = round(float(linja[1]),4)

            kInd = np.where(kRange==k)[0]
            sInd = np.where(sigRange==sigma)[0]

            resultsMat[kInd,sInd] = float(linja[2])

        fid.close()

        return kRange,sigRange,resultsMat


if __name__ == "__main__":
    print "Running..."
