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

from htsint.database import db_connect,ask_upass,fetch_annotations
from htsint.database import Taxon,Gene,Uniprot,GoTerm,GoAnnotation

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

    def testUniprot(self):
        """
        test the accession table
        """

        uniprotQuery = self.session.query(Uniprot).filter_by(uniprot_id='P07663').first()
        self.assertEqual(uniprotQuery.uniprot_entry,"PER_DROME")        
        geneQuery = self.session.query(Gene).filter_by(id=uniprotQuery.gene_id).first()
        self.assertEqual(geneQuery.ncbi_id,'31251')

    def testGoTerm(self):
        """
        test the GoTerm table
        """

        termQuery = self.session.query(GoTerm).filter_by(go_id="GO:0007623").first()
        print 'aspect',termQuery.aspect
        print 'description',termQuery.description
        print 'name',termQuery.name
        
        self.assertEqual(termQuery.aspect,"biological_process")
        self.assertEqual(termQuery.name,"circadian rhythm")
        descLook = re.search("that recurs with a regularity of approximately 24 hours.",termQuery.description)
        self.assertTrue(termQuery.description,descLook)

    def testGoAnnotation(self):
        """
        test the GoTerm table
        """
    
        print("fetching annotations for uniprot id...")
        annotations1 = fetch_annotations(['P07663'],self.session,idType='uniprot')
        self.assertTrue("circadian rhythm" in [a.name for a in annotations1['P07663']])

        print("fetching annotations for ncbi gene id...")
        annotations2 = fetch_annotations(['31251'],self.session,idType='ncbi')
        self.assertTrue("circadian rhythm" in [a.name for a in annotations2['31251']])


### Run the tests
if __name__ == '__main__':
    unittest.main()
