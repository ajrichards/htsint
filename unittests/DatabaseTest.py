#!/usr/bin/env python
"""
Database specific tests
The database must be up and running for these tests to pass
See /src/database/HOWTO
"""

import sys,os,unittest,time,re,time
from sqlalchemy.sql import select
from htsint.database import db_connect,ask_upass,fetch_annotations,fetch_taxa_annotations
from htsint.database import Taxon,Gene,Uniprot,GoTerm,GoAnnotation

## global variables
UPASS = ask_upass()

## test class for the main window function
class DatabaseTest(unittest.TestCase):
    """
    Run a number of tests using taxa id 5476
    """

    def setUp(self):
        """
        connect to the database
        """

        self.session, self.engine = db_connect(upass=UPASS)
        self.conn = self.engine.connect()
        self.testID = '5476'
    
    def testTaxa(self):
        """
        test the taxa table
        """

        query = self.session.query(Taxon).filter_by(ncbi_id=self.testID).all() 
        self.assertEqual(len(query),1)
        query = query[0]
        self.assertTrue(query not in [None])
        self.assertEqual(int(query.ncbi_id),int(self.testID))
        self.assertEqual(query.name,"Candida albicans")
        
    def testGene(self):
        """
        test the gene table
        """

        ## orm
        query1 = self.session.query(Gene).filter_by(ncbi_id='359').first()
        self.assertEqual(query1.symbol,'AQP2')
        self.assertEqual(self.session.query(Taxon).filter_by(id=query1.taxa_id).first().ncbi_id,9606)
        
        ## core
        myDict = {}
        s = select([Gene.symbol]).where(Gene.ncbi_id=='359')
        _result = self.conn.execute(s)
        result = _result.fetchall()
        self.assertEqual(result[0]['symbol'],'AQP2')

    def testUniprot(self):
        """
        test the uniprot table
        """

        uniprotQuery = self.session.query(Uniprot).filter_by(uniprot_ac='P41181').first()
        uniprotGeneQuery = self.session.query(Gene).filter_by(id=uniprotQuery.gene_id).first()
        
        self.assertEqual(uniprotQuery.uniprot_ac,"P41181")
        self.assertEqual(uniprotQuery.uniprot_entry,"AQP2_HUMAN")
        self.assertEqual(uniprotGeneQuery.ncbi_id,"359")
        self.assertEqual(self.session.query(Taxon).filter_by(id=uniprotQuery.taxa_id).first().ncbi_id,9606)
        self.assertEqual(self.session.query(Taxon).filter_by(id=uniprotGeneQuery.taxa_id).first().ncbi_id,9606)

    def testGoTerm(self):
        """
        test the GoTerm table
        """

        termQuery = self.session.query(GoTerm).filter_by(go_id="GO:0007623").first()
        self.assertEqual(termQuery.aspect,"biological_process")
        self.assertEqual(termQuery.name,"circadian rhythm")
        descLook = re.search("that recurs with a regularity of approximately 24 hours.",termQuery.description)
        self.assertTrue(termQuery.description,descLook)

    def testFetchAnnotations(self):
        """
        test the GoTerm and GoAnnotation tables
        """
    
        print("fetching annotations for uniprot id...")
        annotations1 = fetch_annotations(['P41181'],self.engine,idType='uniprot',useIea=False,verbose=True)
        termNames1 = [a[1] for a in annotations1['P41181']]
        self.assertTrue("water transport" in termNames1)

        print("fetching annotations for ncbi gene id...")
        annotations2 = fetch_annotations(['359'],self.engine,idType='ncbi',useIea=False)
        termNames2 = [a[1] for a in annotations2['359']]
        self.assertTrue("water transport" in termNames1)
   
    def testFetchTaxaAnnotations(self):
        """
        test the GoTerm and GoAnnotation tables
        """
    
        print("fetching annotations for taxa")
        geneAnnots,uniprotAnnots = fetch_taxa_annotations([self.testID],self.engine,useIea=False,verbose=True)
        self.assertTrue('GO:0018343' in uniprotAnnots['Q9Y765'])
        
    '''
    def testIdUniqueness(self):
        """
        test that gene ids are unique and we have only one entry for each
        test that the uniprot entries are unique and we have only one for each
        """
        
        taxonId = '9606'
        taxaQuery = self.session.query(Taxon).filter_by(ncbi_id=taxonId).first()
        geneIds = [g.ncbi_id for g in self.session.query(Gene).filter_by(taxa_id=taxaQuery.id).all()]
        self.assertEqual(len(geneIds), len(list(set(geneIds))))

        uniprotIds = [u.uniprot_ac for u in self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()]
        self.assertEqual(len(uniprotIds), len(list(set(uniprotIds))))
    '''    

    '''
    def testAnnotationEquality(self):
        """
        Fetching annotations with uniprot + gene entities should get the same results as with taxa_id
        This test takes some time but should still pass 
        commented out by default
        """

        taxonId = '7091'
        taxaQuery = self.session.query(Taxon).filter_by(ncbi_id=taxonId).first()
        annotations1 = self.session.query(GoAnnotation).filter_by(taxa_id=taxaQuery.id).all()

        geneIds = [g.id for g in self.session.query(Gene).filter_by(taxa_id=taxaQuery.id).all()]
        uniprotIds = [u.id for u in self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()]

        geneAnnotations = self.session.query(GoAnnotation).filter(GoAnnotation.gene_id.in_(geneIds)).all()
        uniprotAnnotations = self.session.query(GoAnnotation).filter(GoAnnotation.uniprot_id.in_(uniprotIds)).all()
        annotations2 = list(set(geneAnnotations).union(set(uniprotAnnotations)))

        self.assertEqual(len(annotations1),len(annotations2))

    '''

### Run the tests
if __name__ == '__main__':
    unittest.main()
