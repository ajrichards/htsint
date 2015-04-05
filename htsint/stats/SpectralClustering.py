#!/usr/bin/env python
"""
A spectral clustering algorithm 

Originally proposed by Andrew Ng et al.

A. Ng, M. Jordan and Y. Weiss 
On spectral clustering: Analysis and an algorithm
In Advances in Neural Information Processing Systems 14, 2001

Modified for self-tuning of sigma and k by Zelnik-Manor et al.

L. Zelnik-manor and P. Perona
Self-tuning spectral clustering
In Advances in Neural Information Processing Systems 17, 2004

"""

__author__ = "Adam Richards"

import sys,os,re,csv
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
    dtype - 'points','distance','similarity' [default]

    """

    def __init__(self,M,dtype='distance'):
        """
        M - matrix 
        dtype - matrix type [points|distance|similarity|]
        """

        ## error check
        dtype = dtype.lower()
        if dtype not in ['points','distance','similarity']:
            raise Exception("Invalid dtype %s"%dtype)

        ## if M is a filepath with lines (i,j,[dist|sim]) then save/load
        if type(M) == type('str') and os.path.exists(M) and re.search(".csv",M):
            distancesPath = M
            matrixPath = re.sub("\.csv","-matrix.npy",distancesPath)
            genesPath = re.sub("\.csv","-genes.npy",distancesPath)
            if not os.path.exists(matrixPath):
                print("...saving matrix in binary format for faster indexing")
                fid = open(distancesPath,'r')
                reader = csv.reader(fid)
                header = reader.next()

                items = set([])
                for linja in reader:
                    items.update(linja[:1])
                fid.close()
                items = np.sort(np.array(list(items)))

                ## create the similarity matrix
                M = np.zeros((items.size,items.size),)
                fid = open(distancesPath,'r')
                reader = csv.reader(fid)
                header = reader.next()
                for linja in reader:
                    i = np.where(items == linja[0])
                    j = np.where(items == linja[1])
                    M[i,j] = float(linja[2])
                    M[j,i] = float(linja[2])

                fid.close()
                np.save(genesPath,items)
                np.save(matrixPath,M)
            else:
                M = np.load(matrixPath)
                items = np.load(genesPath)

        ## error check
        if type(M) != type(np.array([1,2,3])):
            raise Exception("M must be a matrix or a valid file path to a matrix")

        if dtype == 'points':
            _M = M.copy()
            _M = pdist(_M,'sqeuclidean')
            M = squareform(_M)
            dtype == 'distance'

        if dtype == 'similarity':
            M = np.sqrt(M.max()-M)

        ## add eps to avoid divide by zero
        M = M + np.spacing(1)

        self.items = items
        self.M = M

    def run(self,k,sk=None,sigma=None,verbose=False):
        """
        run spectral clustering
        given a number of clusters (k) and bandwidth param (sigma) 
        
        """

        ## create the affinity matrix
        if verbose == True:
            print "\tsimilarity to affinity matrix..."

        k = int(round(k))

        self.k = k
        self.sigma = sigma

        self.A = self.similarity_to_affinity(self.M,sk=sk,sigma=sigma)
        
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
        
    def similarity_to_affinity(self,dMat,sk=7,sigma=None):
        """
        transform a similarity matrix into an affinity matrix

        """

        ## Ng et al. method
        if sk == None:
            if sigma == None:
                raise Exception("If using Ng et al method must specify sigma")
        
            A = np.exp(-1.0 * (self.M**2.0)  /  2.0 * (sigma**2.0))


        ## using method proposed by Zelnik-Manor et al. method
        else:
            A = np.zeros(dMat.shape)

            if dMat == None:
                print "ERROR: distance matrix is None cannot find affinity"
                return None

            ## get sigma k for each row
            sigmaK = np.zeros(dMat.shape[0])
            for i in range(dMat.shape[0]):
                row = dMat[i,:]
                sigmaK[i] = np.sort(row)[sk-1]

            ## there is prob a faster was of doing this
            for i in range(dMat.shape[0]):
                for j in range(dMat.shape[0]):
                    A[i,j] = np.exp((-1.0 * (dMat[i,j]**2.0))  /  (sigmaK[i] * sigmaK[j]))

        
        ## ensure diag is zeros
        A = A - np.diag(np.diag(A))

        return A


    def save(self,labelsPath='sc-labels.csv'):
        """
        save the labels
        """

        ## save the results 
        outFid = open(labelsPath,'wa')
        writer = csv.writer(outFid)
        writer.writerow(['k=%s'%(self.k),'sigma=%s'%(self.sigma)])
        header = ['gene','label']
        writer.writerow(header)

        for i,gene in enumerate(self.items):
            writer.writerow([self.items[i],self.labels[i]])

        outFid.close()
        print("labels saved (sigma=%s,k=%s)."%(self.sigma,self.k))

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
    sc = SpectralCluster(dataScatter,dtype='points')
    sc.run(2,sk=7)

    ax = fig.add_subplot(131)
    for c in np.unique(dataScatterLabels):
        inds = np.where(dataScatterLabels == c)[0]
        ax.plot([dataScatter[inds,0]],[dataScatter[inds,1]],marker='o',color=colors[c],markersize=8.0)
    for c in np.unique(sc.labels):
        inds = np.where(sc.labels == c)[0]
        ax.plot([dataScatter[inds,0]],[dataScatter[inds,1]],marker=markers[c],color='k',markersize=8.0)        

    ax.set_aspect(1./ax.get_data_ratio())
    
    ## circle plot
    sc = SpectralCluster(dataCircle,dtype='points')
    sc.run(2,sk=7)

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
    sc = SpectralCluster(dataLetters,dtype='points')
    sc.run(3,sk=7)

    ax = fig.add_subplot(133)
    for c in np.unique(dataLettersLabels):
        inds = np.where(dataLettersLabels == c)[0]
        ax.plot([dataLetters[inds,0]],[dataLetters[inds,1]],marker='o',color=colors[c],markersize=8.0)
    for c in np.unique(sc.labels):
        inds = np.where(sc.labels == c)[0]
        ax.plot([dataLetters[inds,0]],[dataLetters[inds,1]],marker=markers[c],color='k',markersize=8.0)        

    ax.set_ylim([3,4])
    ax.set_aspect(1./ax.get_data_ratio())

    fig.savefig("spectral-clustering-test.png")



