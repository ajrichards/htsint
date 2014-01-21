#!/usr/bin/env python
"""
A generic template
"""

__author__ = "Adam Richards"

import numpy as np
from scipy.spatial.distance import pdist,cdist,squareform
from scipy.cluster.vq import kmeans2

from htsint.stats import get_silhouette_values

class SpectralCluster(object):
    """
    Class to perform spectral clustering per the algorithm by A. Ng.

    M - similarity matrix
    k - the number of clusters
    sigma - bandwidth parameter for the kernel

    """

    def __init__(self,M,similarity=True):
        

        if similarity == False:
            _M = M.copy()
            _M = pdist(_M,'sqeuclidean')
            M = squareform(_M)
            M = 1.0/(1.0+M)
            
        self.M = M

    def run(self,k,sigma):
        """
        run spectral clustering
        given a number of clusters (k) and bandwidth param (sigma) 
        
        """

        ## create the affinity matrix
        self.A = self.similarity_to_affinity(self.M,sigma)
        
        ## create the diagonal matrix D
        self.D = np.diag(self.A.sum(axis=1)**-0.5)
        
        ## create the L matrix
        _L = np.dot(self.D,self.A)                      # multiply A and D^{-1/2}
        self.L = np.dot(_L,self.D)                      # multiply the above result times D^{-1/2}
     
        # ensure L contains real, finite numbers
        testNan = np.where(np.isnan(self.L) == True)
        testFinite = np.where(np.isfinite(self.L) == False)

        if np.array([len(z) for z in testFinite]).sum() > 0:
            print "WARNING: failed finite check"
        elif np.array([len(z) for z in testNan]).sum() > 0:
            print "WARNING: failed nan check"

        ## find the k largest eigenvectors of L
        eigVals, eigVecs = np.linalg.eig(self.L)

        eigVecs = -1.0 * eigVecs
        sortedEigInds = np.argsort(np.sum(abs(eigVecs),0))
        self.X = eigVecs[:,sortedEigInds[-k:]]

        ## compute normalized matrix Y from X
        n,k = np.shape(self.X)
        self.Y = np.zeros([n,k])
        unitLengths = np.sum(self.X**2,axis=0)**(0.5)

        for col in range(k):
            self.Y[:,col] = self.X[:,col] / unitLengths[col]

        ## cluster the rows of Y using Kmeans
        tries = 0
        iters = 0
        minDistortion = 1e08
        bestClusters = None
        bestKmeanLabels = None
        repeats = 5

        
        bestRepeat = (None,None,-2.0)
        for repeat in range(repeats):
            kmeanResults, kmeanLabels = kmeans2(self.Y,k,minit='points')
            silValues = get_silhouette_values([self.Y],[kmeanLabels],subsample=10000,
                                                      minNumEvents=3,resultsType='raw')
            avgSilVal = silValues['0'].mean()

            if kmeanResults == None or avgSilVal == -2.0:
                continue

            if avgSilVal > bestRepeat[2]:
                bestRepeat = (kmeanResults,kmeanLabels,avgSilVal)

        if bestRepeat[0] == None:
            raise Exception('Kmeans on Y matrix -- did not obtain results')
            return None

        labels = bestRepeat[1]
        avgSilVal = bestRepeat[2]

   
    def similarity_to_affinity(self,dMat,sigma):
        """
        transform a similarity matrix into an affinity matrix
        """

        if dMat == None:
            print "ERROR: distance matrix is None cannot find affinity"
            return None

        A = np.exp(-1.0 * (self.M**2.0)  /  2.0 * (sigma**2.0))

        return A


if __name__ == "__main__":
    print "Running..."

    from htsint.stats import dataScatter,dataScatterLabels
    from htsint.stats import dataCircle,dataCircleLabels
    from htsint.stats import dataLetters,dataLettersLabels
    from htsint.stats import dataNetwork,dataNetworkLabels

    import matplotlib.pyplot as plt

    ## setup
    colors = ['blue','orange','red','green','yellow','magenta','cyan','black']
    fig = plt.figure()
    
    ## scatter plot
    sc = SpectralCluster(dataScatter,similarity=False)
    sc.run(2,0.1)

    ax = fig.add_subplot(121)
    for c in np.unique(dataScatterLabels):
        inds = np.where(dataScatterLabels == c)[0]
        ax.plot([dataScatter[inds,0]],[dataScatter[inds,1]],marker='o',color=colors[c],markersize=8.0)
    ax.set_aspect(1./ax.get_data_ratio())

    
    
    ax = fig.add_subplot(122)
    for c in np.unique(dataCircleLabels):
        inds = np.where(dataCircleLabels == c)[0]
        ax.plot([dataCircle[inds,0]],[dataCircle[inds,1]],marker='o',color=colors[c],markersize=8.0)
    ax.set_aspect(1./ax.get_data_ratio())

    fig.savefig("spectral-clustering-test.png")
