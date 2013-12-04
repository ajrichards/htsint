#!/usr/bin/env python
"""
Blast class specific tests
you must fetch databases before running blast

See the README.md file

"""

import sys,os,unittest,time,re,time
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('agg')

from htsint import Blast,ParseBlast

## test class for the main window function
class BlastTest(unittest.TestCase):
    """
    Run a number of tests using taxa id 7227
    """

    def setUp(self):
        """
        connect to the database
        """
        
        self.queryFile = "opuntia.fasta"
        self.blast = Blast(self.queryFile)

    def testGetQueryFile(self):
        """
        test the function breaks the fasta file in to chunks
        """
        start,stop = 0,2
        newQueryFile = self.blast.get_query_file(".",start,stop)
        queryFileName = os.path.split(self.queryFile)[-1]
        queryFilePath = os.path.join(".",re.sub("\.\w+","",queryFileName,flags=re.IGNORECASE)+"-%s-%s.fasta"%(start,stop))
        self.assertTrue(os.path.exists(queryFilePath))
        os.remove(queryFilePath)

    def testRunBlastX(self):
        """
        test running blastx
        """
     
        ## run the blast
        outFile = "opuntia-0-2.xml"
        targetDB = "swissprot"
        start,stop = 0,2
        self.blast.run_blastx(targetDB,evalue=0.001,start=start,stop=stop)
    
        ## read the blast
        parser = ParseBlast(outFile,resultsDir=".")
        parser.run()


        ## clean up
        for fn in ["opuntia-0-2_1.csv", "opuntia-0-2_1.log","opuntia-0-2.fasta","opuntia-0-2.xml"]:
            self.assertTrue(os.path.exists(fn))
            os.remove(fn)

### Run the tests
if __name__ == '__main__':
    unittest.main()
