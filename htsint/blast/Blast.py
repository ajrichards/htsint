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

An example of blastx usage is:

blastx -query myfile.fa -db nr -out myfile.xml -evalue 0.001 -outfmt 5

"""

import sys,os,time,re,getopt
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline

class Blast(object):
    """
    A generic class to run blast
    """

    def __init__(self,queryFile,BLASTDB="/usr/local/share/htsint"):
        """
        Constructor
            queryFile - is a fasta file of sequences
        """
        
        ## check the query file
        self.queryFile = os.path.realpath(queryFile)
        if os.path.exists(self.queryFile) == False:
            raise Exception("ERROR: Could not locate the fasta file... \n%s"%self.queryFile)

        ## set the environmental variable if specified
        if BLASTDB != None:
            os.environ['BLASTDB'] = BLASTDB

    def get_query_file(self,outDir,start,stop):
        """
        indexing start = 0 and stop = 2 will return a file with
        the first and second sequences in it.

        This means that we start at index start and we stop at 
        the stop index with inclusion.

        """

        queryFileName = os.path.split(self.queryFile)[-1]
        newQueryFile = os.path.join(outDir,re.sub("\.\w+","",queryFileName,flags=re.IGNORECASE)+"-%s-%s.fasta"%(start,stop))
        handleIn = open(self.queryFile, "rU")
        handleOut = open(newQueryFile, "w")
        toWrite = []

        ## get the total number of queries
        

        ## we only include the start if it is 0
        indices = range(start,stop)

        total = -1
        for record in SeqIO.parse(handleIn,"fasta") :
            total += 1
            if total in indices:
                toWrite.append(record)

        SeqIO.write(toWrite,handleOut,'fasta')
        handleIn.close()
        handleOut.close()

        return newQueryFile

    def run_blast(self,targetDB,outDir=".",evalue=0.05,
                   start=None,stop=None,cmd='blastx'):
        """
        targetDB - is the protein database to blast against
        outDir   - a place to put all output files other than cwd
        evalue   - is the blast evalue
        start    - index of the first seq
        stop     - index of the last seq

        """
        
        ## error checking
        outDir = os.path.realpath(outDir)
        if os.path.isdir(outDir) == False:
            raise Exception("Output directory does not exist\n%s"%outDir)
        if start != None and stop == None:
            raise Exception("'start' and 'stop' must be specified together")
        if stop != None and start == None:
            raise Exception("'start' and 'stop' must be specified together")

        outFileName = os.path.split(self.queryFile)[-1]
        outFile = re.sub("\.\w+","",outFileName,flags=re.IGNORECASE)+".xml"

        ## specify the query file and out file
        if start != None or stop != None:
            outFilePath = os.path.join(outDir,re.sub("\.xml","",outFile,flags=re.IGNORECASE)+"-%s-%s.xml"%(start,stop))
            query = self.get_query_file(outDir,start,stop)
        else:
            outFilePath = os.path.join(outDir,outFile)
            query = self.queryFile

        print 'Running blast'
        timeStart = time.time()
        blastx_cline = NcbiblastxCommandline(query=query, db=targetDB, evalue=evalue,
                                             outfmt=5, out=outFilePath,cmd=cmd)
        stdout, stderr = blastx_cline()
        print "Total run time: %s"%time.strftime('%H:%M:%S', time.gmtime(time.time()-timeStart))

        return outFilePath

if __name__ == "__main__":

    ## read in input file 
    if len(sys.argv) < 3:
        print sys.argv[0] + "-f first -l last -q query_file -d database -e evalue -o outdir"
        sys.exit()

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'f:l:q:d:e:o:c:')
    except getopt.GetoptError:
        raise Exception(sys.argv[0] + "-f first -l last -q query_file -d database -e evalue -o outdir -c cmd")
        sys.exit()

    first,last,query,database,evalue,cmd = None,None,None,None,None,None
    for o,a in optlist:
        if o == '-f':
            first = int(a)
        if o == '-l':
            last  = int(a)
        if o == '-q':
            query = a
        if o == '-d':
            database  = a
        if o == '-e':
            evalue  = float(a)
        if o == '-o':
            outdir  = a
        if o == '-c':
            cmd  = a

    blast = Blast(query)
    blast.run_blast(database,outDir=outdir,start=first,stop=last,evalue=evalue,cmd=cmd)
