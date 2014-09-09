#!/usr/bin/env python
"""
Blast class specific tests
you must fetch databases before running blast

See the README.md file

"""

import sys,os,unittest,time,re,shutil
from Bio import SeqIO
from htsint.blast import Blast,ParseBlast,ParallelBlast
from htsint import __basedir__

## test class for the main window function
class BlastTest(unittest.TestCase):
    """
    Run a number of tests using taxa id 7227
    """

    def setUp(self):
        """
        connect to the database
        """

        self.queryFile = os.path.join(os.path.dirname(__file__),"adh.fasta")

    def testGetQueryFile(self):
        """
        test the function breaks the fasta file in to chunks
        """

        self.blast = Blast(self.queryFile)

        ## index 0 to 2 should return 2 results [0,1]
        start,stop = 0,2
        newQueryFile = self.blast.get_query_file(".",start,stop)
        queryFileName = os.path.split(self.queryFile)[-1]
        queryFilePath = os.path.join(".",re.sub("\.\w+","",queryFileName,flags=re.IGNORECASE)+"-%s-%s.fasta"%(start,stop))
        self.assertTrue(os.path.exists(queryFilePath))

        handleIn = open(newQueryFile, "rU")
        total = 0
        for record in SeqIO.parse(handleIn,"fasta") :
            total += 1

        self.assertEqual(total,2)
        handleIn.close()
        os.remove(queryFilePath)

    def testRunBlastX(self):
        """
        test running blastx
        """
     
        ## run the blast
        self.blast = Blast(self.queryFile)
        outFile = "adh-0-2.xml"
        targetDB = "swissprot"
        start,stop = 0,2
        self.blast.run_blastx(targetDB,evalue=0.1,start=start,stop=stop,cmd='blastx')
    
        ## read the blast
        parser = ParseBlast(outFile,outDir=".")
        parser.run()

        ## clean up
        for fn in ["adh-0-2_1.csv", "adh-0-2_1.log","adh-0-2.fasta","adh-0-2.xml"]:
            self.assertTrue(os.path.exists(fn))
            os.remove(fn)
    
    def testParallelBlast(self):
        """
        test that we can create shell scripts to be run in a parallel env
        """

        parBlast = ParallelBlast(self.queryFile,'swissprot') 
        parBlast.evalue = 0.05
        chunks = 3
        parBlast.create_scripts(chunks,"myemail@somewhere.edu")
        #parBlast.submit()
        shutil.rmtree(os.path.join(".","cluster"))

### Run the tests
if __name__ == '__main__':
    unittest.main()
