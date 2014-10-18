#!/usr/bin/env python
"""
Database specific tests
The database must be up and running for these tests to pass
See /src/database/HOWTO
"""

import sys,os,unittest,time,re,time
from sqlalchemy.sql import select
from htsint.database import db_connect,ask_upass,fetch_annotations,fetch_taxa_annotations
from htsint.database import Taxon,Gene,Refseq,Uniprot,GoTerm,GoAnnotation

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
        self.conn = self.engine.connect()
        self.testID = '7227'

    def testTaxa(self):
        """
        test the taxa table
        """

        query = self.session.query(Taxon).filter_by(ncbi_id=self.testID).all() 
        self.assertEqual(len(query),1)
        query = query[0]
        self.assertTrue(query not in [None])
        self.assertEqual(int(query.ncbi_id),int(self.testID))
        self.assertEqual(query.name,"Drosophila melanogaster")
        self.assertTrue("fruit fly" in [query.common_name_1,query.common_name_2,query.common_name_3])
        
    def testGene(self):
        """
        test the gene table
        """

        ## orm
        self.assertTrue(self.session.query(Gene).count() > 4)       
        query1 = self.session.query(Gene).filter_by(ncbi_id='3771877').first()
        self.assertEqual(query1.symbol,'Adh')
        self.assertEqual(self.session.query(Taxon).filter_by(id=query1.taxa_id).first().ncbi_id,7227)

        query2 = self.session.query(Gene).filter_by(ncbi_id='31251').first()
        self.assertEqual(query2.symbol,'per')
        self.assertEqual(self.session.query(Taxon).filter_by(id=query2.taxa_id).first().ncbi_id,7227)
        
        ## core
        myDict = {}
        s = select([Gene.symbol]).where(Gene.ncbi_id=='3771877')
        _result = self.conn.execute(s)
        result = _result.fetchall()
        self.assertEqual(result[0]['symbol'],'Adh')

    def testRefseq(self):
        """
        test the refseq table
        if there are more than one isoform for a protein we should get back more than one result
        """

        ## orm
        refseqQuery = self.session.query(Refseq).filter_by(protein_accession='NP_001027267.1').all()
        actualGeneQuery = self.session.query(Gene).filter_by(ncbi_id='3771877').first()
        self.assertTrue(len(refseqQuery)>0)
        self.assertTrue(refseqQuery[0].gene_id,actualGeneQuery.id)
        
    def testUniprot(self):
        """
        test the uniprot table
        """

        uniprotQuery = self.session.query(Uniprot).filter_by(uniprot_ac='P07663').first()
        uniprotGeneQuery = self.session.query(Gene).filter_by(id=uniprotQuery.gene_id).first()

        self.assertEqual(uniprotQuery.uniprot_ac,"P07663")
        self.assertEqual(uniprotQuery.uniprot_entry,"PER_DROME")
        self.assertEqual(uniprotGeneQuery.ncbi_id,"31251")
        self.assertEqual(self.session.query(Taxon).filter_by(id=uniprotQuery.taxa_id).first().ncbi_id,7227)
        self.assertEqual(self.session.query(Taxon).filter_by(id=uniprotGeneQuery.taxa_id).first().ncbi_id,7227)

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
        annotations1 = fetch_annotations(['P56645'],self.engine,idType='uniprot',useIea=False,verbose=True)
        termNames1 = [a[1] for a in annotations1['P56645']]
        self.assertTrue("regulation of circadian sleep/wake cycle, sleep" in termNames1)

        print("fetching annotations for ncbi gene id...")
        annotations2 = fetch_annotations(['31251'],self.engine,idType='ncbi',useIea=False)
        termNames2 = [a[1] for a in annotations2['31251']]
        self.assertTrue("circadian rhythm" in termNames2)

    def testFetchTaxaAnnotations(self):
        """
        test the GoTerm and GoAnnotation tables
        """
    
        print("fetching annotations for taxa")
        annotations,goTerms = fetch_taxa_annotations(['7091'],self.engine,idType='uniprot',useIea=False,verbose=True)
        termNames = [goTerms[annot] for annot in annotations['Q9NL89']]
        self.assertTrue("regulation of innate immune response" in termNames)

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

    '''
    def testAnnotationSpeed(self):
        """
        test the GoTerm table
        commented out by default
        """
    
        #taxonId = 1297549
        #taxaQuery = self.session.query(Taxon.id,Taxon.ncbi_id).filter_by(ncbi_id=taxonId).first()
        #uniprotQuery = self.session.query(Uniprot.id).filter_by(taxa_id=taxaQuery.id).all()
        #geneIds = [g.id for g in self.session.query(Gene).filter_by(taxa_id=taxaQuery.id).all()]
        #geneIds = [g for g in self.session.query(Gene.id,Gene.ncbi_id).filter_by(taxa_id=taxaQuery.id).all()]
        #print taxaQuery 
        #print geneIds
        #print len(uniprotQuery)

        print 'blah'
        timeStart = time.time()
        print("fetching annotations for uniprot id...")
        annotations1 = fetch_annotations(['P56645','O54943'],self.engine,self.session,idType='uniprot',asTerms=False)
        terms1 = [self.session.query(GoTerm).filter_by(id = a.go_term_id).first().name for a in annotations1['P07663']]
        print terms1
        print("fetch time: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
        #self.assertTrue("circadian rhythm" in terms1)
        
        #print("fetching annotations for ncbi gene id...")
        #annotations2 = fetch_annotations(['31251'],self.session,idType='ncbi',asTerms=False)
        #terms2 = [self.session.query(GoTerm).filter_by(id = a.go_term_id).first().name for a in annotations2['31251']]
        #self.assertTrue("circadian rhythm" in terms2)
    '''


### Run the tests
if __name__ == '__main__':
    unittest.main()
