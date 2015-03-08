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

        self.taxId = '5476'
        
    def testCheckTaxon(self):
        """
        ensure taxon check works
        """

        go = GeneOntology([self.taxId],upass=UPASS)
        go.check_taxon(self.taxId)
    
    def test01GetDicts(self):
        """
        ensure we can create, save and retrieve gene2go and go2gene dictionaries
        """

        termsPath = 'terms.pickle'
        if os.path.exists(termsPath) == True:
            os.remove(termsPath)

        go = GeneOntology([self.taxId],upass=UPASS,idType='ncbi',useIea=True,\
                          aspect='biological_process')
        go.create_dicts(termsPath)
        gene2go, go2gene = go.load_dicts(termsPath)
        print("there are %s genes"%(len(gene2go.keys())))
        print("there are %s terms"%(len(go2gene.keys())))

    
    def test02CreateGoGraph(self):
        """
        ensure we can create, save and retrieve the gograph
        """

        termsPickle = 'terms.pickle'
        graphPickle = 'graph.pickle'
        
        go = GeneOntology(self.taxId,upass=UPASS,idType='ncbi',useIea=True)
        G = go.create_gograph(termsPath=termsPickle,graphPath=graphPickle)

        for picklePath in [termsPickle,graphPickle]:
            if os.path.exists(picklePath):
                os.remove(picklePath)
    
### Run the tests
if __name__ == '__main__':
    unittest.main()
