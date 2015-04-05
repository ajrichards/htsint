.. pipeline example

Gene set analysis
======================

The procedure in this tutorial example to carry out gene set analysis using ``htsint`` is composed of three main steps. 

   #. Gene set generation
   #. Gene set testing
   #. Gene set visualization

This example takes approximately 30 minutes to complete on a modern desktop computer.  The code discussed in this document is available as a script for convenience.

   * :download:`gsa-example.py`

Gene set generation
----------------------------

The basic process involves integrating Gene Ontology [Ashburner00]_ information from one or more taxa to infer functional distances between genes.  These distances are then used to cluster the genes in a specified list.  These genes are any groups of genes that you might want to cluster.  For example, you may use the genes in a large pathway or all genes from an RNA-Seq experiment, as is the case in this example.

Create a term graph
^^^^^^^^^^^^^^^^^^^^^^^^^^

Because genes and their ontology terms will be loaded multiple times fetch the annotations only once then save the dictionaries.

   >>> from htsint import GeneOntology
   >>> go = GeneOntology(["8364","8355"],useIea=False,aspect='molecular_function')
   >>> termsPath = "go-terms.pickle"
   >>> graphPath = "go-graph.pickle"
   >>> go.create_dicts(termsPath)
   >>> gene2go,go2gene = go.load_dicts(termsPath)
   >>> G = go.create_gograph(termsPath=termsPath,graphPath=graphPath)
   >>> print("%s genes have at least one annotation"%(len(gene2go.keys())))
   Term graph for with 15719 nodes successfully created.
   >>> print("Term graph for with %s nodes successfully created."%(len(G.nodes())))
   1291 genes have at least one annotation

The valid GO aspects are: 'biological_process', 'molecular_function' and 'cellular_component'.

Calculate term distances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ideally, this step is carried out in a cluster environment and if you are using `Grid Engine <http://gridscheduler.sourceforge.net>`_ then there are built-in convenience methods.  Whether you are in a high performance environment or on a single machine the initialization is the same. 

   >>> from htsint import TermDistances
   >>> td = TermDistances(termsPath,graphPath)
   >>> termDistancePath = "term-distances.csv"

Using Grid Engine:
"""""""""""""""""""""

   >>> cpus = 60
   >>> td.create_scripts('youremail@somewhere.edu',cpus=cpus)
   >>> td.submit()

Before you submit you can check in the ``htsint-tmp`` directory that was created in the current working directory to ensure the Bash scripts work for your computing environment.  The results are then assembled into a single file.

   >>> from htsint import AssembleDistances
   >>> ad = AssembleDistances(termsPath,resultsPath=termDistancePath)
   >>> ad.run(cpus=cpus)

Using single machine
""""""""""""""""""""""

   >>> td.run_with_multiprocessing(termDistancePath,cpus=16)

This is the most computationally expensive step in the pipeline so for lists with more than a few thousand genes this calculation becomes difficult outside of a cluster environment.  Using 16 cores on a single machine the previous command finished in 00:27:32 (hh:mm:ss).

Calculate gene distances
^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the term-term distances stored in the distance file we can map the gene-gene distances.

   >>> from htsint import GeneDistances
   >>> geneDistancePath = "gene-distances.csv"
   >>> gd = GeneDistances(termsPath,graphPath,termDistancePath,outFile=geneDistancePath)
   >>> gd.run()

Spectral Clustering
^^^^^^^^^^^^^^^^^^^^^^^^^

With the gene-gene distances a number of unsupervised clustering algorithms can be used here.  Because spectral clustering is appropriate for networks we have implemented two algorithms as part of ``htsint``.  There is a bandwidth parameter :math:`\sigma` and a parameter for the number of clusters `k` that need to be given.

Parameter estimation [optional]
"""""""""""""""""""""""""""""""""

   >>> from htsint.stats import SpectralClusterParamSearch, SpectralClusterResults
   >>> scps = SpectralClusterParamSearch(geneDistancePath,dtype='distance')
   >>> scps.run(chunks=15)
   >>> psFigureFile = os.path.join(homeDir,"param-scan-%s.png"%(_aspect))
   >>> scr = SpectralClusterResults(silvalFile,clustersFile)
   >>> scr.plot(figName=psFigureFile,threshMin=5,threshoMax=100)

.. figure:: ../figures/param-scan-mf.png
   :scale: 30%
   :align: center
   :alt: top 75 transcripts
   :figclass: align-center

Ideally, we are looking for values of :math:`\sigma` and `k` that maximize our silhouette value, while at the same time maximize the number of clusters that fall into a reasonable size range.  The size range can be set with the ``threshMin`` and ``threshMax`` arguments.  It helps result interpretation if the specified range can be reasonably investigated through visualization.  The top panel shows the average silhouette value for the clustering results over a grid of possible parameter values. For the same grid the bottom panel illustrates the percentage of total genes that fall into clusters of the desired size.  There is usually a trade-off between high silhouette values and the reasonably sized clusters.  The top three optimal values are marked on the plots.  For this example the parameters are maximized at :math:`k=123` and :math:`\sigma=0.08`.  It is worth noting that strongly associated clusters tend to remain mostly intact over a wide range of parameter values.  In the script version of this example this section the parameter estimation is commented out to minimize compute time.


Run spectral clustering
"""""""""""""""""""""""""""""""""

There are two implementations of spectral clustering available through the SpectralCluster class.  If the argument ``sk`` is ``None`` then the original algorithm proposed by Andrew Ng *et al*. is used [Ng01]_.  Alternatively, a self-tuning version of this algorithm was proposed by Zelnik-Manor and Perona that uses a different :math:`\sigma` around each neighborhood.  The neighborhood size is controlled by the parameter ``sk`` as discussed in the manuscript [Zelnik-Manor04]_.  For smaller networks the self-tuning method gives reasonable results, however for larger networks the grid parameter search seems to provide more biologically intuitive clusters.

   >>> from htsint.stats import SpectralCluster
   >>> k = 123
   >>> sigma = 0.08
   >>> sc = SpectralCluster(geneDistancePath,dtype='distance')
   >>> sc.run(k,sk=None,sigma=sigma,verbose=True)
   >>> sc.save(labelsPath=labelsPath)

Save gene sets
^^^^^^^^^^^^^^^^^^^^

1. Run :doc:`BLAST and create a summarized blast map <blast>`.  To save time in this tutorial we provide an example summary file below.

   * :download:`blast-parsed-summary.csv <blast-parsed-summary.csv>`

   Load the file.

   >>> from htsint.blast import BlastMapper
   >>> bm = BlastMapper()
   >>> bmap = bm.load_summary('blast-parsed-summary.csv',best=False)


   >>> from htsint import GeneSetCollection
   
