#!/usr/bin/env python
"""
A class to represent and draw gene sets
"""

__author__ = "Adam Richards"

from htsint.database import db_connect,Gene,Taxon
from htsint.blast import BlastMapper

class GeneSet(object):
    """
    gene set class
    """


    def __init__(self,geneList):
        """
        Constructor

        geneList

        """

        ## setup db
        self.session, self.engine = db_connect()
        self.conn = self.engine.connect()

        ## declare variables
        self.geneList = list(set(geneList))
        self.geneInfo = self.get_gene_info()

        for gene in self.geneList:
            print gene, self.geneInfo[gene]

    def get_gene_info(self):
        print("...getting gene info")
        geneInfo = {}
        results = self.conn.execute(Gene.__table__.select(Gene.ncbi_id.in_(self.geneList)))
        for row in results:
            geneInfo[str(row.ncbi_id)] = {'symbol': str(row.symbol),
                                          'description': str(row.description),
                                          'taxa_id': self.session.query(Taxon).filter_by(id=row.taxa_id).first().ncbi_id}
        
        return geneInfo

        #geneIds = [g.id for g in self.session.query(Gene).filter_by(taxa_id=taxaQuery.id).all()]     
        #uniprotIds = [u.id for u in self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()]
        #geneAnnotations = self.session.query(GoAnnotation).filter(GoAnnotation.gene_id.in_(geneIds)).all() 
        


if __name__ == "__main__":
    print "Running..."
