.. main file for lpedit documentation

BLAST
======================

BLAST or Basic Local Alignment Search Tool is a widely used tool in computational biology.  A general introduction is beyond the scope of these documents, but more information can be found by visiting the following sites.

   * `NCBI BLAST page <http://blast.ncbi.nlm.nih.gov/Blast.cgi>`_
   * `NCBI BLAST FAQs <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ>`_
   * `BLAST QuickStart <http://www.ncbi.nlm.nih.gov/books/NBK1734/>`_ 

Install and setup
-----------------------------

Local BLAST can be run through ``htsint`` or through the command line. BLAST results help integrate data across taxa through comparative genomics.  For local BLAST functionality install the Suite of `BLAST+ tools <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_.  Under Ubuntu you can simply use:

   .. code-block:: bash

      ~$ sudo apt-get install ncbi-blast+

BLAST
-----------------

If you have successfully :doc:`setup your database <database>` then the `SwissProt database <http://web.expasy.org/docs/swiss-prot_guideline.html>`_ has been downloaded into your data directory.  Any other database or organism that you wish to BLAST against may be downloaded to this directory.  For example, to download the protein sequences of *Anolis carolinensis* and make the appropriate BLAST files use the following.

   .. code-block:: bash
 
      ~$ cd /your/data/directory
      ~$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Anolis_carolinensis/protein/protein.fa.gz
      ~$ gunzip -c protein.fa.gz > Anolis_carolinensis.fasta
      ~$ makeblastdb -in Anolis_carolinensis.fasta -dbtype 'prot' -out Anolis_carolinensis 

To run BLAST from within ``htsint`` use the Blast class.  For convenience ``htsint`` uses the sequence classes in `BioPython <http://biopython.org/wiki/Main_Page>`_ to handle basic operations.  Here is the FASTA file used in this example.

   * :download:`adh.fasta <../data/adh.fasta>`

>>> from htsint.blast import Blast
>>> blast = Blast('adh.fasta')
>>> outFile = "adh-0-2.xml"
>>> targetDB = "uniprot_sprot.fasta.db"
>>> blast.run_blast(targetDB,evalue=0.0001,start=0,stop=2,cmd='blastx')
Running blast
Total run time: 00:00:30

In this example the output file contains only two sequences as constrained by start and stop. Indexes in Python begin with zero.  If you type ``range(start,stop)`` it can help you understand exactly which indexes are implied.  The start and stop options are optional, but they are convenient when moving BLAST into a cluster environment.  If stop is larger than the number of sequences it will still stop at the last one in the file.

This produces the XML output file that was specified and a log file.  To do BLAST the entire file from the command line use you may use:

   .. code-block:: bash

      export BLASTDB='/usr/local/share/htsint/:$BLASTDB'
      blastx -query adh.fasta -db uniprot_sprot.fasta.db -out adh.xml -evalue 0.0001 -num_threads 8 -outfmt 5

Because BLAST results are fairly large by default we parse them into a more reasonable size.

>>> from htsint.blast import ParseBlast
>>> parser = ParseBlast(outFile)
>>> parser.run()
total hits: 2
sequences with at least one match : 2
sequences without any matches: 0
complete.

BLAST in a cluster environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wrappers for the Blast and ParseBlast classes have been created for grid engine cluster environments.

>>> from htsint.blast import ParallelBlast
>>> parBlast = ParallelBlast(queryFile.fasta,'uniprot_sprot.fasta.db')
>>> parBlast.evalue = 0.0001
>>> chunks = 30
>>> parBlast.create_scripts(chunks,"myemail@somewhere.edu")
>>> parBlast.submit()

The individual scripts will be created in a folder called `cluster` and chunks are the number of jobs you wish to break the FASTA file into.  To parse the results use:

>>> from htsint.blast import ParseParallelBlast
>>> parser = ParseParallelBlast(queryFilePath)
>>> chunks = 30
>>> parser.run(chunks)


BlastMapper
---------------------

BlastMapper is a class that htsint uses to interact with BLAST results.  BlastMapper takes as input the parsed results.  A method called ``create_summarized`` will take the parsed results and try to find matching UniProt and NCBI identifiers.  This summary can take time so it is saved as a separate file.


>>> from htsint.blast import BlastMapper

