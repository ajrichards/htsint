#!/usr/bin/env python

import os,sys
import numpy as np
from SilValueGenerator import SilValueGenerator

def get_silhouette_values(matList,matLabelList,subsample=None,minNumEvents=4,resultsType='clusterMeans'):
    '''
    returns a dict of results where files are indexed by a string of the int index
    resultsType -- defines the value of
    example $> silResults = get_silhouette_values([matData1, matData2], [matLabels1, matLabs2])
    silResults 
    '''

    silValues = {}
    silValuesElements = {}
    numFiles = len(matLabelList)
    for expName in range(numFiles):
        silValues[str(expName)] = {}

    ## create subset if data for large data sets 
    subsetExpData = []
    subsetExpLabels = []

    if subsample != None:
        for expInd in range(numFiles):
            expData = matList[expInd]
            expLabels = matLabelList[expInd]
            newIndices = []

            totalInds = 0
            for cluster in np.sort(np.unique(expLabels)):
                clusterInds = np.where(expLabels==cluster)[0]
                totalInds += len(clusterInds)

                if len(clusterInds) > subsample:
                    percentTotal = float(len(clusterInds)) / float(len(expLabels)) 
                    randSelected = np.unique(clusterInds[np.random.randint(0,len(clusterInds),subsample)])
                    newIndices += randSelected.tolist()
                else:
                    newIndices += clusterInds.tolist()

            ## save indices and data
            subsetExpData.append(expData[newIndices,:])
            subsetExpLabels.append(expLabels[newIndices])

    ## calculate the silvalues for each file and the subsampled clusters
    for fileInd in range(numFiles):
            
        if subsample != None:
            fileData = subsetExpData[fileInd]
            fileLabels = subsetExpLabels[fileInd]
        else:
            fileData = matList[fileInd]
            fileLabels = matLabelList[fileInd]

        fileClusters = np.sort(np.unique(fileLabels))    
    
        ## calculate sil values
        svg = SilValueGenerator(fileData,fileLabels)
        silValuesElements[str(fileInd)] = svg.silValues
        
        ## save only sil values for each cluster
        if resultsType == 'clusterMeans':
            for clusterID in fileClusters:
                clusterElementInds = np.where(fileLabels == clusterID)[0]
                if len(clusterElementInds) < minNumEvents:
                    silValues[str(fileInd)][str(clusterID)] = None
                else:
                    clusterSilValue = silValuesElements[str(fileInd)][clusterElementInds].mean()
                    silValues[str(fileInd)][str(clusterID)] = clusterSilValue

                del clusterElementInds
        elif resultsType == 'raw':
            silValues[str(fileInd)] = svg.silValues

    return silValues
