#!/usr/bin/python
"""

"""

import sys
import numpy as np
from htsint.tools import read_matrix,read_de_results,Heatmap

## load differential expression data
deseqIds, deseqColumns, deseqMat = read_de_results('deseq.csv',tool='DESeq')
dfeIds,dfeColumns,dfeMat = read_matrix('deseq-samples.csv',mtype='float')
padjInd = np.where(deseqColumns == 'padj')[0]

## filter out nans 
print deseqMat.shape
padjInd = np.where(deseqColumns == 'padj')[0]
size1 = deseqIds.shape[0]
filter1 = np.where(~np.isnan(deseqMat[:,padjInd]))[0]
deseqIds = deseqIds[filter1]
deseqMat = deseqMat[filter1,:]
mask = np.in1d(dfeIds,deseqIds)
dfeIds = dfeIds[mask]
dfeMat = dfeMat[mask,:]
print("... %s/%s transcripts pass nan filter"%(filter1.size,size1))

## filter for only the most significant transcripts (max 50)
threshold = 0.5
size2 = deseqIds.shape[0]
filter2 = np.where(deseqMat[:,padjInd] <= threshold)[0][:50]
deseqIds = deseqIds[filter2]
deseqMat = deseqMat[filter2,:]
mask = np.in1d(dfeIds,deseqIds)
dfeIds = dfeIds[mask]
dfeMat = dfeMat[mask,:]
print("... %s/%s transcripts pass significance (%s) filter"%(filter2.size,size2,threshold))

## create the heatmap
rowLabels = dfeIds
colLabels = dfeColumns
hm = Heatmap(dfeMat,rowLabels,colLabels)
hm.draw_heatmap(cmap='uy',clabels=True,rlabels=True,rowFont=6)
hm.save("heatmap-demo.png",dpi=200)
