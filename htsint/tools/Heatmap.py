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

def blue_black_yellow():
    cdict = {'red':   ((0.0, 0.0, 0.0),\
                       (0.5, 0.0, 0.1),\
                       (1.0, 1.0, 1.0)),\

             'green': ((0.0, 0.0, 0.0),\
                       (0.5, 0.1, 0.0),\
                       (1.0, 1.0, 1.0)),\

             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def red_black_green():
    cdict = {'red': ((0.0, 0.0, 0.0),\
                     (0.5, 0.0, 0.1),\
                     (1.0, 1.0, 1.0)),\

             'green': ((0.0, 0.0, 1.0),\
                       (0.5, 0.1, 0.0),\
                       (1.0, 0.0, 0.0)),

             'blue': ((0.0, 0.0, 0.0),\
                      (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def red_black_blue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),

             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

class Heatmap(object):
    """
    produce a clustering heatmap

    linkage metrics = 'average', 'single', 'centroid', 'complete'

    """

    def __init__(self,mat,rowLabels,colLabels,width=6,height=7,title=None,dpi=300,
                 hpad=0.14,vpad=0.05,fontSize=10,fontName='sans-serif'):
        """ 
        Constructor
        """

        self.indx = {}
        self.z = {}
        self.fontSize = fontSize
        self.fontName = fontName
        self.mat = mat

        if type(rowLabels) == type([]):
            self.rowLabels = np.array(rowLabels)
        else:
            self.rowLabels = rowLabels
        
        if type(colLabels) == type([]):
            self.colLabels = np.array(colLabels)
        else:
            self.colLabels = colLabels

        ## initialize the axis (l,b,w,h)
        self.plt = plt
        self.fig = self.plt.figure(figsize=(width,height),dpi=dpi)

        ## set relative size of subplots (sum to 1.0)
        base1 = 0.15
        base2 = 0.85

        ## top dendogram
        self.ax1 = self.fig.add_axes([base1,vpad+base2,base2-hpad,base1-vpad])
        self.ax1.set_yticks([])
        self.ax1.set_xticks([])
        self.ax1.set_frame_on(False)

        ## side dendogram
        self.ax2 = self.fig.add_axes([0,vpad,base1,base2])
        self.ax2.set_yticks([])
        self.ax2.set_xticks([])
        self.ax2.set_frame_on(False)

        ## heatmap
        self.ax3 = self.fig.add_axes([base1,vpad,base2-hpad,base2])
        self.ax3.set_yticks([])
        self.ax3.set_xticks([])
        self.ax3.set_frame_on(False)

        ## colorbar
        hsize = base1           # same as side dendogram width
        vsize = base1-vpad      # same as top dendogram height
        left = 0.1 * hsize
        bottom = vpad+base2 + (0.4 * vsize)
        width = 0.80 * hsize
        height = 0.5 * vsize
        self.ax4 = self.fig.add_axes([left,bottom,width,height],frame_on=False)
        self.ax4.set_yticks([])
        self.ax4.set_xticks([])
        self.ax4.set_frame_on(False)

        self.title = None

        ## cluster the matrix
        self.mat = mat
        self.cluster(0)
        self.cluster(1)

    def cluster(self,dim,labels=None):
        """
        use hierarchical clustering to group the data
        dim = 0 are the rows
        dim = 1 are the columns
        """
        
        if dim == 0:
            print("clustering the rows...")
            x = self.mat
            orientation='right'
            self.plt.sca(self.ax2)
            ax = self.ax2
        if dim == 1:
            print("clustering the columns...")
            x = self.mat.transpose()
            orientation='top'
            self.plt.sca(self.ax1)
            ax = self.ax1

        distMatrix = pdist(x)
        distMatrix = squareform(distMatrix)
        linkageMatrix = linkage(distMatrix,method='complete')

        z = dendrogram(linkageMatrix,orientation=orientation,\
                       no_labels=True,color_threshold=1.0,\
                       link_color_func=lambda k: 'k')

        indx = z['leaves']
        self.indx[str(dim)] = indx
        self.z[str(dim)] = z

    def draw_heatmap(self,cmap='uy',clabels=True,rlabels=False,rowFont=None,colFont=None):
        """
        draw the heatmap portion of the plot
        cmap can be a custom instance of a cmap
        or 'yu' for yellow-black-blue'
        or 'rg' for red-black-green
        """
        
        if cmap == 'rg':
            cmap = red_black_green()
        elif cmap == 'uy':
            cmap = blue_black_yellow()
        else:
            cmap = cmap

        if rowFont == None:
            rowFont = self.fontSize-2
        if colFont == None:
            colFont = self.fontSize

        self.plt.sca(self.ax3)
        ax = self.ax3

        n,m = self.mat.shape
        
        if not self.indx.has_key('0') or not self.indx.has_key('1'):
            raise Exception("cluster before plotting heatmap")

        ## setup event handler
        def mat_picker(x, mouseevent):
            """ 
            find the row and column
            """

            if mouseevent.xdata == None or mouseevent.ydata == None:
                return False, dict()

            colCoord = np.floor(mouseevent.xdata+.5)
            rowCoord = np.floor(mouseevent.ydata+.5)
            
            
            if len(self.colLabels) > 0:
                col = self.colLabels[self.indx['1']][colCoord]
            else:
               col = np.arange(m)[self.indx['1']][colCoord]

            if len(self.rowLabels) > 0:
                row = self.rowLabels[self.indx['0']][rowCoord]
            else:
                row = np.arange(n)[self.indx['0']][rowCoord]

            print("row: %s (%s), col: %s (%s)"%(col,int(colCoord),row,int(rowCoord)))
            return False, dict()

        ## reorder matrix
        matReordered = self.mat[self.indx['0'],:]
        matReordered = matReordered[:,self.indx['1']]
        hmap = ax.imshow(matReordered, interpolation='nearest',aspect='auto',cmap=cmap,picker=mat_picker,origin='lower')

        ## handle axes
        if clabels and len(self.colLabels) > 0:
            ax.set_xticks(range(m))
            ax.set_xticklabels(self.colLabels[self.indx['1']],fontsize=colFont,fontname=self.fontName,rotation='vertical')
        
        if rlabels and len(self.rowLabels) > 0:
                
            if n > 200:
                print ("WARNING: too many rows to visualize row labels")
                rowFont = None

            if rowFont:
                ax.yaxis.set_ticks_position('right')
                ax.set_yticks(range(n))
                ax.set_yticklabels(self.rowLabels[self.indx['0']],fontsize=rowFont,fontname=self.fontName)

        ## colorbar
        self.plt.sca(self.ax4)
        ax = self.ax4

        val = np.ceil(np.abs(self.mat).max())
        norm = mpl.colors.Normalize(vmin=-1*val, vmax=val)
        cb1 = mpl.colorbar.ColorbarBase(ax,cmap=cmap,
                                        ticks=[int(round(i)) for i in np.linspace(-1*val,val,5)],
                                        norm=norm,
                                        orientation='horizontal')

        for t in ax.get_xticklabels():
            t.set_fontsize(self.fontSize)
            t.set_fontname(self.fontName)
        for t in ax.get_yticklabels():
            t.set_fontsize(self.fontSize)
            t.set_fontname(self.fontName)

        #print "min: %s,max: %s,mean: %s"%(round(mat.min(),2),round(mat.max(),2),round(mat.mean(),2))

    def save(self,fileName,dpi=300):
        """
        save the current plot to file
        """

        self.plt.savefig(fileName,dpi=dpi)

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
    mat[:,3:] = mat[:,3:] -2.0

    hm = Heatmap(colLabels=np.array(["A","B","C","D","E","F"]),\
                 rowLabels= np.array(["r"+str(i) for i in range(n*2)]))
    hm.cluster(mat,0)
    hm.cluster(mat,1)
    hm.draw_heatmap(mat,cmap='uy',clabels = True, rlabels=True)

    ## error checking
    hm.show()
