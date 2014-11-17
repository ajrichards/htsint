#!/usr/bin/env python
"""
produce a clustering heatmap

http://docs.scipy.org/doc/scipy-0.14.0/reference/cluster.hierarchy.html

"""

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram


__author__ = "Adam Richards"

class Heatmap(object):
    """
    produce a clustering heatmap
    """

    def __init__(self,width=6,height=6,title=None,dpi=400):
        """ 
        Constructor
        """

        ## initialize the axis (l,b,w,h)
        self.plt = plt
        self.fig = self.plt.figure(figsize=(width,height),dpi=dpi)

        self.buff1 = 0.0001
        self.buff2 = 0.0001
        
        ## top dendogram
        self.ax1 = self.fig.add_axes([0.3,self.buff2+.7,0.7,0.3])
        self.ax1.set_yticks([])
        self.ax1.set_xticks([])

        ## side dendogram
        self.ax2 = self.fig.add_axes([self.buff1,self.buff2,0.3,0.7])
        self.ax2.set_yticks([])
        self.ax2.set_xticks([])
        
        ## heatmap
        self.ax3 = self.fig.add_axes([0.3,self.buff2,0.7,0.7])
        self.ax3.set_yticks([])
        self.ax3.set_xticks([])

        ## colorbar
        self.ax4 = self.fig.add_axes([0.018,.815,.263,0.1],frame_on=False)
        self.ax4.set_yticks([])
        self.ax4.set_xticks([])

        self.title = None

    def cluster(self,mat,dim):
        """
        use hierarchical clustering to group the data
        dim = 0 are the rows
        dim = 1 are the columns
        """
        
        if dim == 0:
            print("clustering the rows...")
            x = mat
            orientation='right'
            self.plt.sca(self.ax2)
            ax = self.ax2
        if dim == 1:
            print("clustering the columns...")
            x = mat.transpose()
            orientation='top'
            self.plt.sca(self.ax1)
            ax = self.ax1

        distMatrix = pdist(x)
        distMatrix = squareform(distMatrix)
        linkageMatrix = linkage(distMatrix,method='complete')
        
        dendrogram(linkageMatrix,orientation=orientation,
                   color_threshold=0.3,leaf_font_size=6)

    def draw_heatmap(self,mat,cmap=None):
        """
        draw the heatmap portion of the plot
        """
        
        if cmap == None:
            cmap = self.plt.cm.PuBuGn

        self.plt.sca(self.ax3)
        ax = self.ax3

        n,m = mat.shape
        hmap = ax.imshow(mat, interpolation='nearest',aspect='auto',cmap=cmap)

        if self.title != None:
            ax.set_title(self.title)

        #ax.set_xticks(range(m))
        #ax.set_yticks(range(n))
        #ax.set_xticklabels(["c"+str(i) for i in range(m)])
        #ax.set_yticklabels(["r"+str(i) for i in range(n)])

        ## colorbar
        self.plt.sca(self.ax4)
        ax = self.ax4
        norm = mpl.colors.Normalize(vmin=-3, vmax=3)
        
        cb1 = mpl.colorbar.ColorbarBase(ax,cmap=cmap,
                                        ticks=[-2, 0, 2],
                                        norm=norm,
                                        orientation='horizontal')

        print "min: %s,max: %s,mean: %s"%(round(mat.min(),2),round(mat.max(),2),round(mat.mean(),2))

    def save(self,fileName):
        """
        save the current plot to file
        """

        self.plt.savefig("example_heatmap.pdf")

    def show(self):
        """
        display the current plot
        """

        self.plt.show()

if __name__ == "__main__":
    print "Running..."


    n = 10
    m = 6

    mat = np.vstack((np.random.normal(0,1,(n,m)),np.random.normal(3,1,(n,m))))
    mat[:,3:] = mat[:,3:] + 1.0

    hm = Heatmap()
    hm.cluster(mat,0)
    hm.cluster(mat,1)
    hm.draw_heatmap(mat)
    hm.show()
    
