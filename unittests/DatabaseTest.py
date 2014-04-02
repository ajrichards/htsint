#!/usr/bin/env python
"""
Database specific tests
The database must be up and running for these tests to pass
See /src/database/HOWTO
"""

import sys,os,unittest,time,re,time
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('agg')

from htsint.database import db_connect,ask_upass
from htsint.database import Taxon,Gene,Accession,GoTerm,GoAnnotation

## global variables
UPASS = ask_upass()


## test class for the main window function
class DatabaseTest(unittest.TestCase):
    """
    Run a number of tests using taxa id 7227
    """

    def setUp(self):
        """
        connect to the database
        """

        self.session, self.engine = db_connect(upass=UPASS)
        self.testID = '7227'

    def testTaxa(self):
        """
        test the taxa table
        """

        query = self.session.query(Taxon).filter_by(ncbi_id=self.testID).first() 
        self.assertTrue(query not in [None])
        self.assertEqual(int(query.ncbi_id),int(self.testID))
        self.assertEqual(query.name,"Drosophila melanogaster")
        self.assertTrue("fruit fly" in [query.common_name_1,query.common_name_2,query.common_name_3])

    def testGene(self):
        """
        test the gene table
        """

        self.assertTrue(self.session.query(Gene).count() > 4)       
        query = self.session.query(Gene).filter_by(ncbi_id='3771877').first() 
        self.assertEqual(query.symbol,'Adh')
        
    def testAccession(self):
        """
        test the accession table
        """

        geneQuery = self.session.query(Gene).filter_by(ncbi_id='32006').first()
        accessionQuery = self.session.query(Accession).filter_by(gene_id=geneQuery.id).all()
        self.assertTrue('161077725' in [aq.protein_gi for aq in accessionQuery])

    def testGoTerm(self):
        """
        test the GoTerm table
        """

        termQuery = self.session.query(GoTerm).filter_by(go_id="GO:0008150").first()
        self.assertEqual(termQuery.aspect,"Process")
        self.assertEqual(termQuery.description,"biological_process")

    def testGoAnnotation(self):
        """
        test the GoTerm table
        """
    
        geneQuery = self.session.query(Gene).filter_by(ncbi_id='3771877').first() 
        annotationQuery = self.session.query(GoAnnotation).filter_by(gene_id=geneQuery.id).all()
        annotations = [aq.go_term_id for aq in annotationQuery]
        terms = [self.session.query(GoTerm).filter_by(id = a).first().go_id for a in annotations]
        self.assertTrue("GO:0004022" in terms)

### Run the tests
if __name__ == '__main__':
    unittest.main()
