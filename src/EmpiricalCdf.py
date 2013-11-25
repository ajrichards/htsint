#!/usr/bin/env python
from __future__ import division
import re
import numpy as np


###########################################################

class EmpiricalCdf:

    def __init__(self,data,numBins=10):

        if type(data) == type([]):
            data = np.array(data)

        self.data = np.sort(data)

        ## variables
        self.numBins = numBins

        ## def get range and frequencies
        if re.search('int',str(type(self.data[0]))):
            self.dRange,self.freqs = self.get_discrete_bins(self.data)
            self.isDiscrete = True
        elif re.search('float',str(type(self.data[0]))):
            self.dRange,self.freqs = self.get_continuous_bins(self.data)
            self.isDiscrete = False
        else:
            print 'ERROR" EmpiricalCDF bad element type', type(self.data[0])
            return None

        ## get the cmf and pmf
        self.pmf, self.cmf = self.get_pmf_and_cmf(self.freqs)

        ## specify percentiles
        totalObservations = float(len(self.data))
        self.percentiles = np.arange(totalObservations + 1)[1:] / totalObservations

    def get_discrete_bins(self,data):
        '''
        --get_bins--
        given a numpy array or a list fn returns
            dRange - a range of all the potential values of the data
            dFreqs - corresponding to dRange all the observed frequencies
                 for the input data 
        '''

        dRange = np.arange(data.min(),data.max()+1)
        dFreqs = np.array(map(lambda i: np.where(data==dRange[i])[0].size,range(dRange.size)))

        return dRange, dFreqs

    def get_continuous_bins(self,data):
        hist,binEdges = np.histogram(data,self.numBins)
        return binEdges, hist

    def get_pmf_and_cmf(self,freqList):
        '''
        --get_pmf_and_cmf--
        given a numpy array of bins (see get_bins)
        return an empirical pmf and a empirical cmf 
        '''

        pmf,cmf = [],[]

        total = freqList.sum()
        for f in range(len(freqList)):
            pmf.append(float(freqList[f]) / total)
            cmf.append(sum(pmf[:f + 1]))

        return pmf,cmf

    def get_percentile(self,x):
        '''
        function to return the percentile associated with a given int x
        '''

        minX = self.dRange.min()
        maxX = self.dRange.max()
        
        if x < minX:
            return 0.0
        elif x > maxX:
            return 1.0
        else:
            if self.isDiscrete == True:
                indX = np.where(self.dRange==x)[0][0]
            else:
                indX = 0
                for i in range(len(self.dRange)):
                    if x >= self.dRange[i]:
                        indX = i
                        if indX >= len(self.cmf):
                            indX = len(self.cmf) - 1
                        
            return self.cmf[indX]

    def get_value(self,x):
        '''
        returns a value given a percentile
        '''

        minX = self.dRange.min()
        maxX = self.dRange.max()

        if x <= 0:
            return 0.0
        elif x < minX:
            indX = 0

        if x >= 1.0:
            return maxX

        indX = 0
        for i in range(len(self.percentiles)):
            if x >= self.percentiles[i]:
                indX = i

        #print 'debug', len(self.data), len(self.percentiles), indX, float(indX) / float(len(self.percentiles))

        return self.data[indX]

### run the test 
if __name__ == '__main__':
    '''
    can be verified with the ecdf function in R
    '''

    ## discrete test
    a = [1,1,1,2,2,2,3,3,5,5]
    ecdf = EmpiricalCDF(a) 

    if ecdf.get_percentile(0) != 0.0:
        print 'ERROR: 0 value'
    if ecdf.get_percentile(1) != 0.3:
        print 'ERROR: 1 value'
    if ecdf.get_percentile(2) != 0.6:
        print 'ERROR: 2 value'
    if ecdf.get_percentile(3) != 0.8:
        print 'ERROR: 3 value'
    if ecdf.get_percentile(4) != 0.8:
        print 'ERROR: 4 value'
    if ecdf.get_percentile(5) != 1.0:
        print 'ERROR: 5 value'
    if ecdf.get_percentile(100) != 1.0:
        print 'ERROR: 100 value'


    if ecdf.get_value(0.0) != 0:
        print 'ERROR: 0 value'
    if ecdf.get_value(0.3) != 1:
        print 'ERROR: 0.3 value'
    if ecdf.get_value(0.6) != 2:
        print 'ERROR: 0.6 value'
    if ecdf.get_value(0.7) != 3:
        print 'ERROR: 0.7 value'
    if ecdf.get_value(0.8) != 3:
        print 'ERROR: 0.8 value'
    if ecdf.get_value(1.0) != 5:
        print 'ERROR: 1.9 value'
    if ecdf.get_value(100) != 5:
        print 'ERROR: 100 value',

    print '.........'

    ## continuous percentile test
    x = [0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.5,0.5]
    eCDF = EmpiricalCDF(x)
    
    if eCDF.get_percentile(0.0) != 0.0:
        print 'ERROR: 0 value'
    if eCDF.get_percentile(0.1) != 0.3:
        print 'ERROR: 0.1 value'
    if eCDF.get_percentile(0.2) != 0.6:
        print 'ERROR: 0.2 value'
    if eCDF.get_percentile(0.3) != 0.8:
        print 'ERROR: 0.3 value'
    if eCDF.get_percentile(0.4) != 0.8:
        print 'ERROR: 0.4 value'
    if eCDF.get_percentile(0.5) != 1.0:
        print 'ERROR: 0.5 value'
    if eCDF.get_percentile(100) != 1.0:
        print 'ERROR: 100 value'

    ## continuous value test
    if eCDF.get_value(0.0) != 0.0:
        print 'ERROR: 0 value'
    if eCDF.get_value(0.3) != 0.1:
        print 'ERROR: 0.1 value'
    if eCDF.get_value(0.6) != 0.2:
        print 'ERROR: 0.2 value'
    if eCDF.get_value(0.7) != 0.3:
        print 'ERROR: 0.3 value'
    if eCDF.get_value(0.8) != 0.3:
        print 'ERROR: 0.4 value',eCDF.get_value(0.8)
    if eCDF.get_value(1.0) != 0.5:
        print 'ERROR: 0.5 value'
    if eCDF.get_value(100) != 0.5:
        print 'ERROR: 100 value'

    print '0.3', eCDF.get_value(0.3)


    if eCDF.get_value(1.0) != 0.5:
        print 'ERROR: 1.0 value'


    print 'finished test.'
