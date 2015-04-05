#!/usr/bin/env python
"""
A generic template
"""

import os,csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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

    def plot(self,threshMax=100,threshMin=5,fontSize=10,fontName='sans-serif',cmap=plt.cm.PuOr,figName='param-scan.png'):
        """
        create a heatmap plot
        top panel are the sil values
        bottom panel denotes cluster sizes with respect to a specified range 
        """
        
        clustersMat = np.zeros((self.kRange.size,self.sigRange.size),)

        for k in self.kRange:
            for sigma in self.sigRange:
                clusters = self.clusters[str(k)][str(sigma)]
                tooSmall = np.where(clusters < threshMin)[0]
                tooLarge = np.where(clusters > threshMax)[0]
                tooSmallGenes = np.array([clusters[ts] for ts in tooSmall]).sum()
                tooLargeGenes = np.array([clusters[tl] for tl in tooLarge]).sum()
                percentAccepted = ((clusters.sum() - tooSmallGenes - tooLargeGenes) / clusters.sum())
                kInd = np.where(self.kRange==k)[0]
                sInd = np.where(self.sigRange==sigma)[0]
                clustersMat[kInd,sInd] = percentAccepted

        ## get best
        combined = clustersMat + self.silvals
        cols = np.argsort(combined.max(axis=0))[::-1][:3]
        rows = np.argsort(combined.max(axis=1))[::-1][:3]
        print("The maximimum values are:")
        print("best k: %s"%self.kRange[rows[0]])
        print("best sigma: %s"%self.sigRange[cols[0]])  

        ## create the figure
        fig = plt.figure(figsize=(7,6))
        ax1 = plt.subplot2grid((2, 5), (0, 0),colspan=4)
        ax2 = plt.subplot2grid((2, 5), (1, 0),colspan=4)
        ax4 = plt.subplot2grid((2, 5), (0, 4),rowspan=2)    

        ## sil value panel
        ax1.plot(cols,rows,color='k',marker='x',markersize=5,markeredgewidth=4,linestyle='None')
        p1 = ax1.imshow(self.silvals, interpolation='nearest',vmin=-1.0,vmax=1.0, origin='lower',aspect='auto',cmap=cmap)
        ax1.set_xticks(range(self.sigRange.shape[0]))
        ax1.set_yticks(range(self.kRange.shape[0]))
        ax1.set_xticklabels([round(i,2) for i in self.sigRange],rotation=45,fontsize=fontSize,fontname=fontName)
        ax1.set_yticklabels([int(round(i)) for i in self.kRange],fontsize=fontSize,fontname=fontName)
        ax1.set_title("Silhouette values",fontsize=fontSize+2,fontname=fontName)
        ax1.set_ylabel(r"$k$",fontsize=fontSize+1,fontname=fontName)
        ax1.set_xlabel(r"$\sigma$",fontsize=fontSize+1,fontname=fontName)
        #ax1.yaxis.set_major_locator(MaxNLocator(5))

        ## cluster size panel
        ax2.plot(cols,rows,color='k',marker='x',markersize=5,markeredgewidth=4,linestyle='None')
        p2 = ax2.imshow(clustersMat, interpolation='nearest',vmin=-1.0,vmax=1.0, origin='lower',aspect='auto',cmap=cmap)
        ax2.set_xticks(range(self.sigRange.shape[0]))                                                        
        ax2.set_yticks(range(self.kRange.shape[0]))
        ax2.set_xticklabels([round(i,2) for i in self.sigRange],rotation=45,fontsize=fontSize,fontname=fontName)
        ax2.set_yticklabels([int(round(i)) for i in self.kRange],fontsize=fontSize,fontname=fontName)
        ax2.set_title(r"Cluster size $\geq " + str(threshMin) + "$ and $\leq " + str(threshMax) + "$ (%)",fontsize=fontSize+2,fontname=fontName)
        ax2.set_ylabel(r"$k$",fontsize=fontSize+1,fontname=fontName)
        ax2.set_xlabel(r"$\sigma$",fontsize=fontSize+1,fontname=fontName)
        #ax2.yaxis.set_major_locator(MaxNLocator(5)) 

        ## add text
        plt.figtext(0.07,0.92,"A",weight='bold')
        plt.figtext(0.07,0.42,"B",weight='bold')

        ## colorbar
        norm = mpl.colors.Normalize(vmin=-1.0, vmax=1.0)
        cb1 = mpl.colorbar.ColorbarBase(ax4,cmap=cmap,norm=norm,orientation='vertical')
        fig.subplots_adjust(wspace=0.45,hspace=0.5)
        plt.savefig(figName,dpi=400) 

if __name__ == "__main__":
    print "Running..."
