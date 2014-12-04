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


    def testSummarize(self):
        """
        test the summarize function
        """

        self.bm.load_parsed(self.parsedFile)

        print 'cool'

### Run the tests
if __name__ == '__main__':
    unittest.main()
