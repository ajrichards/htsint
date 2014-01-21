#!/usr/bin/env python

import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform

class SilValueGenerator():
    
    def __init__(self,mat,labels,dMetric='euclidean'):
        '''
        mat - is an n x d matrix of observations and features
        labels = is a vector of lenth n with a single assignment for each observation
        dMetric = is a given distance metric

        '''
        
        ## error checking
        matRows,matCols = np.shape(mat)
        if matRows != len(labels):
            print "INPUT ERROR: SilValueGenerator must have matching matrix and labels"

        self.mat = mat
        self.labels = np.array([l for l in labels])
        self.dMetric = dMetric

        ## get values 
        self.silValues = self._get_silhouette_values()


    def _get_silhouette_values(self):

        if self.dMetric != 'euclidean':
            print "ERROR: SilValueGenerator bad distance metric specified"
            return None

        ## get dims
        if len(self.mat.shape) == 1:
            d = 1
            n = self.mat.size
        else:
            n,d = self.mat.shape

        distInter = np.zeros((n,),)
        distIntra = np.zeros((n,),)
        silhouetteVals = None
        uniqueLabels = np.sort(np.unique(self.labels))
    
        ## get distances between cluster centers
        centroids = np.zeros((len(uniqueLabels),d))
        for ul in range(len(uniqueLabels)):
            ulab = uniqueLabels[ul]
            indicesK = np.where(self.labels==ulab)[0]
            centroids[ul,:] = self.mat[indicesK,:].mean(axis=0)

        centroidDMatrix = squareform(pdist(centroids,'euclidean'))
        for ul in range(len(uniqueLabels)):
            centroidDMatrix[ul,ul] = np.inf

        ## calculate distances 
        for ul in range(len(uniqueLabels)):
            ulab = uniqueLabels[ul]

            ## within values
            indicesK = np.where(self.labels==ulab)[0]
            elementsK = self.mat[indicesK,:]
            clustDists = pdist(elementsK,'euclidean')
            clustDists = squareform(clustDists)
            #means = np.zeros(len(indicesK,),)
            #for i in range(len(indicesK)):
            #    rowIndices = list(set(range(len(indicesK))).difference(set([i])))
            #    means[i] = clustDists[rowIndices,:].mean()

            means = clustDists.mean(axis=0)
            distIntra[indicesK] = means
            
            ## get closest cluster
            closeClust = uniqueLabels[np.where(centroidDMatrix[ul,:] == centroidDMatrix[ul,:].min())[0][0]]
            
            ## between values
            indicesClose = np.where(self.labels==closeClust)[0]
            elementsClose = self.mat[indicesClose,:]
            clustDists = cdist(elementsClose, elementsK,'euclidean')
            means = clustDists.mean(axis=0)
            distInter[indicesK] = means

        silhouetteVals = (distInter - distIntra) / np.maximum(distIntra, distInter)

        return silhouetteVals
