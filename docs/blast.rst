.. main file for lpedit documentation

BLAST
======================

BLAST or Basic Local Alignment Search Tool is a widely used tool in computational biology.  A general introduction is beyond the scope of these documents, but more information can be found by visiting the following sites.

   * `NCBI BLAST page <http://blast.ncbi.nlm.nih.gov/Blast.cgi>`_
   * `NCBI BLAST FAQs <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ>`_
   * `BLAST QuickStart <http://www.ncbi.nlm.nih.gov/books/NBK1734/>`_ 

Install and setup
-----------------------------

Local BLAST can be run through ``htsint`` or through the command line. BLAST results help integrate data across taxa through comparative genomics.  For local BLAST functionality, install the Suite of `BLAST+ tools <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_.  Under Ubuntu you can simply use:

   .. code-block:: bash

      ~$ sudo apt-get install ncbi-blast+

BLAST
-----------------

If you have successfully setup the ``htsint`` :doc:`database <database>` then the `SwissProt <http://web.expasy.org/docs/swiss-prot_guideline.html>`_ database has been downloaded into your data directory.  Any other database or organism that you wish to BLAST against may be downloaded to this directory.  For example, to download the protein sequences of *Anolis carolinensis* and make the appropriate BLAST files use the following.

   .. code-block:: bash
 
      ~$ cd /your/data/directory
      ~$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Anolis_carolinensis/protein/protein.fa.gz
      ~$ gunzip -c protein.fa.gz > Anolis_carolinensis.fasta
      ~$ makeblastdb -in Anolis_carolinensis.fasta -dbtype 'prot' -out Anolis_carolinensis 

To run BLAST from within ``htsint`` use the Blast class.  For convenience ``htsint`` uses the sequence classes in `BioPython <http://biopython.org/wiki/Main_Page>`_ to handle basic operations.  Here is the FASTA file used in this example.

   * :download:`xtrop.fasta <../data/xtrop.fasta>`

>>> from htsint.blast import Blast
>>> blast = Blast('xtrop.fasta')
>>> targetDB = "uniprot_sprot.fasta.db"
>>> blastFile = blast.run_blast(targetDB,evalue=0.0001,start=0,stop=2,cmd='blastx')
Running blast
Total run time: 00:00:07
>>> print(blastFile)
./xtrop-0-2.xml

In this example the output file contains only two sequences as constrained by start and stop. Indexing in Python begin with zero.  If you type ``range(start,stop)`` it can help you understand exactly which indexes are implied.  The start and stop options are optional, but they are convenient when moving BLAST into a cluster environment.  If stop is larger than the number of sequences it will still stop at the last one in the file.  The outfile name is automatically created with numbers indicating the first and last sequence of the input FASTA file.

This produces the XML output file that was specified and a log file.  To do BLAST the entire file from the command line use you may use:

   .. code-block:: bash

      ~$ export BLASTDB='/usr/local/share/htsint/:$BLASTDB'
      ~$ blastx -query xtrop.fasta -db uniprot_sprot.fasta.db -out xtrop.xml -evalue 0.0001 -num_threads 8 -outfmt 5

Because BLAST results are fairly large by default we parse them into a more reasonable size.

>>> from htsint.blast import ParseBlast
>>> parser = ParseBlast(blastFile)
>>> parsedFile = parser.run()
total hits: 2
sequences with at least one match : 2
sequences without any matches: 0
complete.
>>> print(parsedFile)
./xtrop-0-2_parsed.csv

BLAST in a cluster environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wrappers for the Blast and ParseBlast classes have been created for grid engine cluster environments.

>>> from htsint.blast import ParallelBlast
>>> parBlast = ParallelBlast(queryFile.fasta,'uniprot_sprot.fasta.db')
>>> parBlast.evalue = 0.0001
>>> chunks = 30
>>> parBlast.create_scripts(chunks,"myemail@somewhere.edu")
>>> parBlast.submit()

The individual scripts will be created in a folder called ``cluster`` and chunks are the number of jobs you wish to break the FASTA file into.  To parse the results use:

>>> from htsint.blast import ParseParallelBlast
>>> parser = ParseParallelBlast(queryFilePath)
>>> chunks = 30
>>> parser.run(chunks)

This will create a parsed and concatenated file in the ``cluster`` directory.

BlastMapper
---------------------

BlastMapper is a class that ``htsint`` uses to interact with BLAST results.  BlastMapper takes as input the parsed results.  A method called ``create_summarized`` will take the parsed results and try to find matching UniProt and NCBI identifiers from the local database.  This summary can take time so it is saved as a separate file.

>>> from htsint.blast import BlastMapper
>>> bm = BlastMapper()
>>> summaryFile = bm.create_summarized(parsedFile)
>>> print(summaryFile)
./xtrop-0-2_parsed_summary.csv

BLAST results will most commonly be used in the summarized form so there is some flexibility in how the results may be loaded.  The object of interest is the dictionary that gets returned.  The keys will always correspond to the query identifiers.

>>> bmap = bm.load_summary(summaryFile,best=True,evalue=0.0001)
>>> print(bmap.keys())
['Xetro.K04500.1', 'Xetro.K04761.1']
>>> print(bmap['Xetro.K04500.1'])
('CJ088_HUMAN', '80007', 'Homo sapiens', '9606', 3.528e-35)

Because the argument ``best`` was set to True only the best matching result for each query is returned.  A result is a tuple of the form

   .. code-block:: none

      (database_id, ncbi_gene_id, species_name, species_id, e-value)

If ``best`` is set to False then all hits passing the threshold e-value (default=0.0001) are returned.  The returned object changes from a tuple to a list of tuples.

>>> bmap = bm.load_summary(summaryFile,best=False,evalue=0.0001)
>>> for hit in bmap['Xetro.K04500.1']:
        print(hit)
('CJ088_HUMAN', '80007', 'Homo sapiens', '9606', 3.528e-35)
('CJ088_MOUSE', '68277', 'Mus musculus', '10090', 1.11258e-31)
('CJ088_RAT', '309029', 'Rattus norvegicus', '10116', 9.44503e-30)
('YJ016_HUMAN', '-', 'Homo sapiens', '9606', 7.09757e-09)

We can also restrict the results by taxa.

>>> bmap = bm.load_summary(summaryFile,best=False,evalue=0.0001,taxaList=['10090','9606'])
>>> for hit in bmap['Xetro.K04500.1']:
        print(hit)
('CJ088_HUMAN', '80007', 'Homo sapiens', '9606', 3.528e-35)
('CJ088_MOUSE', '68277', 'Mus musculus', '10090', 1.11258e-31)
('YJ016_HUMAN', '-', 'Homo sapiens', '9606', 7.09757e-09)
