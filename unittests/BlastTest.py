#!/usr/bin/env python
"""
Blast class specific tests
you must fetch databases before running blast

~$ cd /src/blast/
~$ python FetchBlastDBs.py

"""

import sys,os,unittest,time,re,time
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('agg')

from htsint import Blast

## test class for the main window function
class BlastTest(unittest.TestCase):
    """
    Run a number of tests using taxa id 7227
    """

    def setUp(self):
        """
        connect to the database
        """

        targetFile = "opuntia.fasta"
        self.blast = Blast(targetFile)

    def testGeneLists(self):
        """
        make sure gene lists are properly created via both methods
        """
        
        self.blast.run

        pass

### Run the tests
if __name__ == '__main__':
    unittest.main()
