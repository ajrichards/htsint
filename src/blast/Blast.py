#!/usr/bin/env python

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

    pass
    """
    ## specify the isotig file path
    fastaFilePath = os.path.join(os.getenv("HOME"),'research','mobigen','raw_data',"HP10IZV0-12.adaptorcleaner.454Isotigs.fna")

    if os.path.exists(fastaFilePath) == False:
        print "ERROR: Could not locate the fasta file... exiting"
        print fastaFilePath
        sys.exit()

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
