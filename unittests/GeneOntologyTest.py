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

        self.taxID = 7227
        self.geneList = ['30970','30971','30972','30973','30975','30976','30977','30978','30979','30980',
                         '30981','30982','30983','30984','30985','30986','30988','30990','30991','30994']

    def testGeneLists(self):
        """
        make sure gene lists are properly created via both methods
        """
        
        ## via taxid method
        go = GeneOntology(taxID=self.taxID,upass=UPASS)
        self.assertTrue(len(go.geneList) > 1e04)
        self.assertEqual(go.taxID,self.taxID)

        ## via a gene list (mixed taxa mode)
        go = GeneOntology(geneList=self.geneList,upass=UPASS)
        self.assertEqual(len(go.geneList),len(self.geneList))
        self.assertEqual(go.taxID,None)

    def testCheckTaxon(self):
        """
        ensure taxon check works
        """

        go = GeneOntology(taxID=self.taxID,upass=UPASS)
        go.check_taxon(self.taxID)
    
    def testGetTerms(self):
        """
        fetch term-gene relations from db
        """

        go = GeneOntology(geneList=self.geneList,upass=UPASS)
        gene2go = go.get_terms()
        self.assertTrue(len(gene2go.keys()) > 10)

    '''
    def testGetDicts(self):
        """
        ensure we can create, save and retrieve gene2go and go2gene dictionaries
        """

        dictsPickle = 'foo.pickle'
        if os.path.exists(dictsPickle) == True:
            os.remove(dictsPickle)

        go = GeneOntology(geneList=self.geneList,upass=UPASS)
        gene2go1,go2gene1 = go.get_dicts(filePath=dictsPickle)
        gene2go2,go2gene2 = go.get_dicts(filePath=dictsPickle)

        self.assertEqual(len(gene2go1.keys()),len(gene2go2.keys()))
        self.assertEqual(len(go2gene1.keys()),len(go2gene2.keys()))

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

        go = GeneOntology(geneList=self.geneList,upass=UPASS)
        G = go.create_gograph(termsPath=termsPickle,graphPath=graphPickle)
        
        print 'nodes', len(G.nodes())

        for picklePath in [termsPickle,graphPickle]:
            if os.path.exists(picklePath):
                os.remove(picklePath)
    '''

### Run the tests
if __name__ == '__main__':
    unittest.main()
