#!/usr/bin/python
"""
run after assembleDistances.py
 
example:

 ~$ python runSpectralClustering.py -n x-iea -a bp

"""

__author__ = "Adam Richards"

import sys,os,csv,time,gc,re
from multiprocessing import Pool, cpu_count
import numpy as np
from htsint.stats import SpectralCluster


def mp_worker((k,sigma,distancePath,dtype)):
    """
    run spectral clustering
    """

    sc = SpectralCluster(distancePath,dtype=dtype)
    sc.run(k,sk=None,sigma=sigma,verbose=True)
    
    return sc

class SpectralClusterParamSearch(object):
    """
    Class to perform spectral clustering parameter search.

    a distance or similarity matrix path is required 
    where each row in the file is (header is true)

        i,j,[dist|similarity]

    """

    def __init__(self,distancePath,dtype='distance',aspect='biological_process'):
        """
        distancePath - path to distance matrix
        dype = distance [default] or similarity
        aspect = biological_process, molecular_function, or cellular_component
        """

        ## input 
        if os.path.exists(distancePath) == False:
            raise Exception("cannot find distances file %s\nexiting..."%(distancePath))

        if dtype not in ['distance','similarity']:
            raise Exception("Invalid distant type (dtype) specified")

        ## call an instance of SpectralClustering to ensure labels and matrix files are saved
        self.dtype = dtype
        self.aspect = aspect
        self.distancePath = distancePath
        sc = SpectralCluster(distancePath,dtype=dtype)
        matrixPath = re.sub("\.csv","-matrix.npy",distancePath)
        genesPath = re.sub("\.csv","-genes.npy",distancePath)
        self.M = np.load(matrixPath)
        self.items = np.load(genesPath)

        ## output
        self.resultsPath1 = re.sub("\.csv","-scparams-sv.csv",distancePath)
        self.resultsPath2 = re.sub("\.csv","-scparams-cl.csv",distancePath)
      
    def run(self,kRange=None,sigmaRange=None,chunks=None):
        """
        chunks - log(size(M)) (can be reduced if there are memory issues)
        kRange - the range of k to search
        sigmaRange - the range of sigma to search
        """

        ## run spectral clustering parameter search
        totalCores = cpu_count()
        totalCores = totalCores - 1

        ## specify the ranges
        kRange = np.array([int(round(i)) for i in np.linspace(20,500,15)])

        
        ## different sigma ranges are appropriate for different GO aspects
        if sigmaRange:
            pass
        elif self.aspect == 'biological_process':
            sigmaRange = np.linspace(0.01,1.0,15)
        elif self.aspect == 'molecular_function':
            sigmaRange = np.linspace(1.0,2.0,15)
        elif self.aspect == 'cellular_component':
            sigmaRange = np.linspace(0.05,1.0,15)
        else:
            raise Exception("invalid aspect provided")

        ## prepare outfiles
        outFid1 = open(self.resultsPath1,'wa')
        self.writer1 = csv.writer(outFid1)
        header1 = ['k','sigma','silvalue']
        self.writer1.writerow(header1)
        
        outFid2 = open(self.resultsPath2,'wa')
        self.writer2 = csv.writer(outFid2)
        header2 = ['k','sigma']+range(kRange.max())
        self.writer2.writerow(header2)

        ## limit each iteration to keep memory usage down 
        if chunks:
            pass
        else:
            chunks = int(round((np.log(self.M.shape[0]))))
        print("chunks = %s"%chunks)

        toRun = []
        for k in kRange:
            toRun += [(k,sigma,self.distancePath,self.dtype) for sigma in sigmaRange]

        stopPoints = np.arange(0,len(toRun),chunks)
        if stopPoints[-1] < len(toRun):
            stopPoints = np.hstack([stopPoints[1:],np.array([len(toRun)])])

        begin = 0
        for i,chunk in enumerate(range(stopPoints.size)):
            stop = stopPoints[chunk]
            print('...running %s-%s/%s'%(begin,stop,len(toRun)))
            self.run_sc(toRun,begin,stop)
            begin = stop

        print "complete."
        outFid1.close()
        outFid2.close()

    ## determin all clusters larger than a cutoff and cluster them again
    def get_cluster_sizes(self,sc):
        """
        return cluster sizes from a SpectralClustering object that has results
        """

        clusterSizes = []
        clusters = sc.labels
        for k in np.arange(clusters.max()):
            clusterSizes.append(len(np.where(clusters==k)[0]))  
        
        return clusterSizes

    def run_sc(self,toRun,begin,stop):
        """
        use multiprocessing to run spectral clustering over a range of values
        """

        po = Pool(processes=cpu_count()-1)
        _results = po.map_async(mp_worker,toRun[begin:stop])
        results =  _results.get()
    
        ## write the cluster sizes to file
        for r,sc in enumerate(results):
            k,sigma,dpath,dtype = toRun[begin:stop][r]
            clusterSizes = self.get_cluster_sizes(sc)      
            self.writer1.writerow([k,sigma] + [round(sc.avgSilValue,4)])
            self.writer2.writerow([k,sigma] + clusterSizes)

        po.close()
        po.join()
