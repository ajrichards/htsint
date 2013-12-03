#!/usr/bin/env python
"""
For more information on blast cmd line tools:

    http://www.ncbi.nlm.nih.gov/books/NBK52640/

For more information ofn blast databases

    ftp://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html

Some useful information

    blastp  - protein query against a protein database
    blastn  - nucleotide query against a nucleotide database
    blastx  - nucleotide query, dynamically translated in all six frames, against a protein database
    tblastn - protein query against a nucleotide database dynamically translated in all six frames
    tblastx - nucleotide query, dynamically translated in all six frames
              against a nucleotide database similarly translated

    nr      -  non-redundant protein sequence database with entries from GenPept, Swissprot, PIR, 
               PDF, PDB, and NCBI RefSeq
    nt      -  nucleotide sequence database, with entries
               from all traditional divisions of GenBank, 
               EMBL, and DDBJ excluding bulk divisions 
               (gss, sts, pat, est, htg divisions) and wgs 
                entries. Not non-redundant.


"""


import sys,os,time
import Bio
from Bio import SeqIO
print "\n...."
print "Using Biopython version %s"%Bio.__version__
from Bio.Blast.Applications import NcbiblastxCommandline

class Blast(object):
    """
    A generic class to run blast


    """

    def __init__(self,targetFile):
        """
        Constructor
         
            targetFile - is a fasta file of sequences
        
        """
        
        self.targetFile = targetFile

        if os.path.exists(self.targetFile) == False:
            print "ERROR: Could not locate the fasta file... exiting"
            print fastaFilePath
            sys.exit()


    def run(self,targetDB,outXML=None):
        """
        targetDB - is the database to blast against
        
        """

        print 'Running blast'





    """
    ## specify the isotig file path
    fastaFilePath = os.path.join(os.getenv("HOME"),'research','mobigen','raw_data',"HP10IZV0-12.adaptorcleaner.454Isotigs.fna")


    ## count total sequences
    total = 0
    handle = open(fastaFilePath, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        total+=1

    handle.close()
    print "Total sequences in FASTA file:", total

    timeStart = time.time()
    blastx_cline = NcbiblastxCommandline(query=fastaFilePath, db="uniprot_sprot.db", evalue=0.001,
                                         outfmt=5, out="brassicae-pooled-sprot.xml",cmd='blastx')
    stdout, stderr = blastx_cline()
    print "Total run time: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))
    """
