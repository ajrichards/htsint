.. differential expression example

Differential expression analysis
======================================

In principal, all that is needed to carry out a transcript-level
differential expression analysis is a count matrix where rows
correspond to transcripts and columns are the samples from the
experiment.  There are a number of programs and possible pipelines
that you can use to obtain a count matrix, but a review of the
different methods is beyond the scope of this documentation.  In this
example we start with this count matrix and make a comparison between
two unpaired groups.

The count matrix in this example was created with the following sequence of tools.

   1. *de novo* assembly using Trinity [Grabherr11]_
   2. Read abundances estimated with Sailfish [Patro14]_

Transcript level differential expression analysis
----------------------------------------------------- 

In this example the `R <http://www.r-project.org>`_ package `DESeq2 <http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ is used to compare two frog groups (endurant vs non-endurant) that have distinct phenotypes [Love14]_.  Sailfish estimates counts at the individual isoform level.  It is important that the count matrix is not transformed or normalized when used as input into DESeq2.

   * :download:`raw-counts.csv <../raw-counts.csv>`
   * :download:`run-deseq.R <../run-deseq.R>`

To run the example you use:

   .. code-block:: bash

      ~$ Rscript run-deseq.R raw-counts.csv

For more details on the functions used within the Rscript see the `DESeq2 documentation <http://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf>`_.  This will produce two files.  


  1. A csv file with the transcript level differential expression results (``deseq.csv``)
  2. A csv file with counts that have been transformed with a regularized logarithm (``deseq-samples.csv``)


The second file was created using the function ``rlog`` provided with DESeq2.  From the manual, we see that the function transforms the count data to the log2 scale "in a way that in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size". DESeq2 does this by fitting a model.  To estimate differential expression (first file) let :math:`K_{ij}` be the count for transcript :math:`i` and sample :math:`j`.  DESeq2 uses a GLM to models these counts with a negative binomial distribution.

   .. math::

      K_{ij} &\sim \textrm{NB}(\mu_{ij},\alpha_{i})\\
      \mu_{ij} &= s_{j} q_{ij}\\
      \textrm{log}_{2}(q_{ij}) &= x_{j} \beta_{i}

where :math:`\mu_{ij}` is the fitted mean and :math:`\alpha_{i}` is the gene specific dispersion parameter.   The fitted mean is composed of a sample specific size factor (:math:`s_{j}`) and :math:`q_{ij}` which corresponds to the "expected true concentration of fragments for sample j.  The coefficients :math:`\beta_{i}` give the log2 fold changes for gene :math:`i` for each column of the GLM model matrix :math:`X`.

Creating a Heatmap
-------------------------

The code in this section is available as a script.

   * :download:`deseq-example.py <deseq-example.py>`

1. Load these two files into Python.

   >>> from htsint.tools import read_matrix,read_de_results,Heatmap
   >>> deseqIds, deseqColumns, deseqMat = read_de_results('deseq.csv',tool='DESeq')
   >>> dfeIds,dfeColumns,dfeMat = read_matrix('deseq_samples.csv',mtype='float')
   >>> padjInd = np.where(deseqColumns == 'padj')[0]

2. Filter out NA values

   >>> import numpy as np  
   >>> padjInd = np.where(deseqColumns == 'padj')[0]
   >>> size1 = deseqIds.shape[0]
   >>> filter1 = np.where(~np.isnan(deseqMat[:,padjInd]))[0]
   >>> deseqIds = deseqIds[filter1]
   >>> deseqMat = deseqMat[filter1,:]
   >>> mask = np.in1d(dfeIds,deseqIds)
   >>> dfeIds = dfeIds[mask]
   >>> dfeMat = dfeMat[mask,:]
   >>> print("... %s/%s transcripts pass nan filter"%(filter1.size,size1))

3. Filter for only the most significant transcripts (max 50)

   >>> threshold = 0.5
   >>> size2 = deseqIds.shape[0]
   >>> filter2 = np.where(deseqMat[:,padjInd] <= threshold)[0][:50]
   >>> deseqIds = deseqIds[filter2]
   >>> deseqMat = deseqMat[filter2,:]
   >>> mask = np.in1d(dfeIds,deseqIds)
   >>> dfeIds = dfeIds[mask]
   >>> dfeMat = dfeMat[mask,:]

4. Draw a heatmap of transformed count data
 
   >>> rowLabels = dfeIds
   >>> colLabels = dfeColumns
   >>> hm = Heatmap(dfeMat,rowLabels,colLabels)
   >>> hm.draw_heatmap(cmap='uy',clabels=True,rlabels=True,rowFont=6)
   >>> hm.save("heatmap_demo.png",dpi=200)

.. figure:: ../figures/heatmap-demo.png
   :scale: 30%
   :align: center
   :alt: top 75 transcripts
   :figclass: align-center


.. note:: The ``rowLabels`` can be changed for more meaningful labels by using :doc:`BLAST and BlastMapper <blast>`
