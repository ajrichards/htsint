#!/usr/bin/env python
"""
Blast class specific tests
you must fetch databases before running blast

See the README.md file

"""

import sys,os,unittest,time,re,shutil
from Bio import SeqIO
from htsint.blast import BlastMapper
from htsint import __basedir__

## test class for the main window function
class BlastTest(unittest.TestCase):
    """
    Run a number of tests using taxa id
    """

    def setUp(self):
        """
        connect to the database
        """

        self.parsedFile = os.path.join(os.path.dirname(__file__),"blast-parsed.csv")
        self.bm = BlastMapper()
        
    def test01Summarize(self):
        """
        test the summarize function
        """
    
        summaryFile = re.sub("\.csv","",self.parsedFile)+"_summary.csv"
        if os.path.exists(summaryFile):
            os.remove(summaryFile)
                          
        self.bm.create_summarized(self.parsedFile)
        self.assertTrue(os.path.exists(summaryFile))

    def test02Something(self):
        """
        read in the results summary
        """

        summaryFile = re.sub("\.csv","",self.parsedFile)+"_summary.csv"
        bmap = self.bm.load_summary(summaryFile,taxaList=["10090"])
        self.assertEqual(bmap['GG11117|c2_g1_i1'][0],'INT1_MOUSE')
        
        bmap = self.bm.load_summary(summaryFile,taxaList=["10090"],trinityIsoform=True)
        self.assertEqual(bmap['GG11117|c2_g1'][0],'INT1_MOUSE')

        bmap = self.bm.load_summary(summaryFile,taxaList=["10090"],trinityIsoform=True,best=False)
        self.assertEqual(bmap['GG11117|c2_g1'][0][0],'INT1_MOUSE')
        self.assertEqual(bmap['GG11117|c2_g1'][0][2],0.0)
       
### Run the tests
if __name__ == '__main__':
    unittest.main()
