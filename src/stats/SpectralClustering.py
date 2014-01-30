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
    similarity - default is a similarity matrix

    """

    def __init__(self,M,similarity=True):

        if similarity == False:
            _M = M.copy()
            _M = pdist(_M,'sqeuclidean')
            M = squareform(_M)
            M = 1.0/(1.0+M)
            
        self.M = M

    def run(self,k,sigma=0.1,verbose=False):
        """
        run spectral clustering
        given a number of clusters (k) and bandwidth param (sigma) 
        
        """

        ## create the affinity matrix
        if verbose == True:
            print "\tsimilarity to affinity matrix..."

        self.A = self.similarity_to_affinity(self.M,sigma)
        
        ## create the diagonal matrix D
        if verbose == True:
            print "\tcreating diagonal matrix..."

        self.D = np.diag(self.A.sum(axis=1)**-0.5)
        
        ## create the L matrix
        if verbose == True:
            print "\tcreating laplacian matrix..."

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
        if verbose == True:
            print "\tfinding eigenvalues and eigenvectors..."

        eigVals, eigVecs = np.linalg.eig(self.L)
        eigVecs = -1.0 * eigVecs
        sortedEigInds = np.argsort(np.sum(abs(eigVecs),0))
        self.X = eigVecs[:,sortedEigInds[-k:]]

        ## compute normalized matrix Y from X
        if verbose == True:
            print "\tnormalizing Y matrix..."

        n,k = np.shape(self.X)
        self.Y = np.zeros([n,k])
        unitLengths = np.sum(self.X**2,axis=0)**(0.5)

        for col in range(k):
            self.Y[:,col] = self.X[:,col] / unitLengths[col]

        ## cluster the rows of Y using Kmeans
        if verbose == True:
            print "\tkmeans clustering to get labels..."

        tries = 0
        iters = 0
        minDistortion = 1e08
        bestClusters = None
        bestKmeanLabels = None
        repeats = 10

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

        self.labels = bestRepeat[1]
        self.avgSilValue = bestRepeat[2]

        return self.avgSilValue
        
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
    colors = ['cyan','orange','red','green','yellow','magenta','blue','black']
    markers = ["x",".","*","^"]
    fig = plt.figure()
    
    ## scatter plot
    sc = SpectralCluster(dataScatter,similarity=False)
    sc.run(2,sigma=10.0)

    ax = fig.add_subplot(131)
    for c in np.unique(dataScatterLabels):
        inds = np.where(dataScatterLabels == c)[0]
        ax.plot([dataScatter[inds,0]],[dataScatter[inds,1]],marker='o',color=colors[c],markersize=8.0)
    for c in np.unique(sc.labels):
        inds = np.where(sc.labels == c)[0]
        ax.plot([dataScatter[inds,0]],[dataScatter[inds,1]],marker=markers[c],color='k',markersize=8.0)        

    ax.set_aspect(1./ax.get_data_ratio())    
    
    ## circle plot
    sc = SpectralCluster(dataCircle,similarity=False)
    sc.run(2,sigma=0.1)

    ax = fig.add_subplot(132)
    for c in np.unique(dataCircleLabels):
        inds = np.where(dataCircleLabels == c)[0]
        ax.plot([dataCircle[inds,0]],[dataCircle[inds,1]],marker='o',color=colors[c],markersize=8.0)
    for c in np.unique(sc.labels):
        inds = np.where(sc.labels == c)[0]
        ax.plot([dataCircle[inds,0]],[dataCircle[inds,1]],marker=markers[c],color='k',markersize=8.0)        

    ax.set_ylim([2.5,7.5])
    ax.set_aspect(1./ax.get_data_ratio())

    ## letters plot
    sc = SpectralCluster(dataLetters,similarity=False)
    sc.run(3,sigma=5)

    ax = fig.add_subplot(133)
    for c in np.unique(dataLettersLabels):
        inds = np.where(dataLettersLabels == c)[0]
        ax.plot([dataLetters[inds,0]],[dataLetters[inds,1]],marker='o',color=colors[c],markersize=8.0)
    for c in np.unique(sc.labels):
        inds = np.where(sc.labels == c)[0]
        ax.plot([dataLetters[inds,0]],[dataLetters[inds,1]],marker=markers[c],color='k',markersize=8.0)        

    ax.set_ylim([3,4])
    ax.set_aspect(1./ax.get_data_ratio())

    print "letters avg sil value", sc.avgSilValue

    fig.savefig("spectral-clustering-test.png")
