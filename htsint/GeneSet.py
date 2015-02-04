#!/usr/bin/env python
"""
A class to represent and draw gene sets
"""

__author__ = "Adam Richards"

from htsint.database import db_connect,Gene,Taxon,GoTerm
from htsint.blast import BlastMapper

class GeneSet(object):
    """
    gene set class
    """

    def __init__(self,geneList,gene2go):
        """
        Constructor

        geneList - list
        gene2go - dictionary

        """

        ## setup db
        self.session, self.engine = db_connect()
        self.conn = self.engine.connect()

        ## declare variables
        self.geneList = list(set(geneList))
        self.geneInfo = self.get_gene_info()
        self.gene2go = gene2go

    def get_gene_info(self):
        print("...getting gene info")
        geneInfo = {}
        results = self.conn.execute(Gene.__table__.select(Gene.ncbi_id.in_(self.geneList)))
        for row in results:
            geneInfo[str(row.ncbi_id)] = {'symbol': str(row.symbol),
                                          'description': str(row.description),
                                          'taxa': self.session.query(Taxon).filter_by(id=row.taxa_id).first().ncbi_id}        
        return geneInfo

    def get_go_terms(self,gene):
        """
        return a human friendly summary of the go terms
        """

        if self.gene2go.has_key(gene):
            terms = [self.session.query(GoTerm).filter(GoTerm.go_id == key).first().name + " (%s)"%key for key in self.gene2go[gene]]
            return ";".join(terms)
        else:
            return 'None'


if __name__ == "__main__":
    print "Running..."
