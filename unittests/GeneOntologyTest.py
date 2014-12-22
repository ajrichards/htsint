#!/usr/bin/env python
"""
GeneOntology class specific tests
The database must be up and running for these tests to pass
See /htsint/database/HOWTO
"""

import sys,os,unittest,time,re,time
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('agg')

from htsint.database import ask_upass
from htsint import GeneOntology

## global variables
UPASS = ask_upass()

## test class for the main window function
class GeneOntologyTest(unittest.TestCase):
    """
    Run a number of tests using taxa id 7227
    """

    def setUp(self):
        """
        simple setup
        """

        self.taxId = 511145
        
    def testCheckTaxon(self):
        """
        ensure taxon check works
        """

        go = GeneOntology([self.taxId],upass=UPASS)
        go.check_taxon(self.taxId)
    
    def testGetDicts(self):
        """
        ensure we can create, save and retrieve gene2go and go2gene dictionaries
        """

        dictsPickle = 'foo.pickle'
        if os.path.exists(dictsPickle) == True:
            os.remove(dictsPickle)

        go = GeneOntology([self.taxId],upass=UPASS,idType='ncbi',useIea=False,\
                          aspect='biological_process')
        gene2go,go2gene = go.get_dicts(termsPath=dictsPickle)

        #self.assertEqual(len(gene2go1.keys()),len(gene2go2.keys()))
        #self.assertEqual(len(go2gene1.keys()),len(go2gene2.keys()))

        if os.path.exists(dictsPickle) == True:
            os.remove(dictsPickle)

    def testCreateGoGraph(self):
        """
        ensure we can create, save and retrieve the gograph
        """

        termsPickle = 'foo1.pickle'
        graphPickle = 'foo2.pickle'
        for picklePath in [termsPickle,graphPickle]:
            if os.path.exists(picklePath) == True:
                os.remove(picklePath)

        go = GeneOntology(self.taxId,upass=UPASS,idType='ncbi')
        G = go.create_gograph(termsPath=termsPickle,graphPath=graphPickle)
        #print 'nodes', len(G.nodes())

        for picklePath in [termsPickle,graphPickle]:
            if os.path.exists(picklePath):
                os.remove(picklePath)

### Run the tests
if __name__ == '__main__':
    unittest.main()
