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

    def __init__(self,silvalFile,clustersFile):
        """
        Constructor
        """

        ## error checking
        for filePath in[silvalFile,clustersFile]:
            if not os.path.exists(filePath):
                raise Exception("could not find file: %s"%filePath)

        clusters = self.load_clusters_file(clustersFile)


    def load_clusters_file(self,clusterFile):
        """
        load the clusters file

        k,sigma,clustid1,clustid2,....clustid_max

        """
        
        fid = open(clusterFile,'r')
        reader = csv.reader(fid)
        header = reader.next()

        for linja in reader:
            print len(linja)

        print len(header)
        #test = reader.next()
        #
        #for i in range(len(test)):
        #    print header[i], test[i]


        print(header)

        fid.close()

if __name__ == "__main__":
    print "Running..."
