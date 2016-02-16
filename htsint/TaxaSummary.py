#!/usr/bin/env python
"""
A class that prints and generates summary information about a species list
with respect to the information contained in the database
"""

__author__ = "Adam Richards"

import csv,re,sys
from htsint.database import db_connect,Taxon,Gene,Uniprot,GoAnnotation,fetch_annotations,get_evidence_codes

class TaxaSummary(object):
    """
    Summarize the information in a taxa list
    """

    def __init__(self,taxaInput,delimiter="\t",upass=''):
        """
        The input is a path to a tab delimited file
        or simply a python list
        The delimiter can also be specified
        The first column is the ncbi taxa id
        All additional columns are optional and user specified 
        """

        self.taxaList = []
        self.session,self.engine = db_connect(upass=upass)
        if type(taxaInput) == type([]):
            for taxon in taxaInput:
                if re.search("\D",taxon):
                    raise Exception("Bad taxa in taxaList (%s)"%(taxon))
                    continue
                self.taxaList.append(taxon)
        else:
            self.taxaList = list(set(self.read_file(taxaInput)))

        print("There are %s unique taxa in the list"%(len(self.taxaList)))

    def print_names_summary(self):
        """
        print to the screen in REstructuredText
        """

        columns = [10,60,40]
        row = "+"
        head = "+"
        for col in columns:
            row += "-"*col+"-+"
            head += "="*col+"=+"

        def show_contents(columns,items,withTrailing=True):
            if len(columns) != len(items):
                raise Exception("Dimension mismatch")

            toPrint = "| "

            for i,item in enumerate(items):
                item = str(item)

                if len(item) >= columns[i]:
                    raise Exception("col %s not large enough min = %s"%(i,len(item)+2))

                toPrint += item+" "*(columns[i]-len(item))+"| "

            print(toPrint[:-1])
            
            if withTrailing:
                print(row)

        print("\nTaxa Summary\n---------------")
        print(row)
        items = ['Taxon ID','Scientific Name','Common Name']
        show_contents(columns,items,withTrailing=False)
        print(head)

        for taxonId in self.taxaList:            
            query = self.session.query(Taxon).filter_by(ncbi_id=taxonId).all()
            if len(query) > 1:
                raise Exception("Taxon id '%s' return multiple values from database"%taxonId)
            if query == None:
                raise Exception("Invalid taxon '%s' provided in list"%taxonId)
            
            show_contents(columns,[query[0].ncbi_id, query[0].name, query[0].common_name_1])

    def print_annotation_summary(self,minCoverage=0.0,useIea=False):
        """
        print a summary of the annotation information
        """
        
        columns = [10,16,18,14,10,10]
        row = "+"
        head = "+"
        for col in columns:
            row += "-"*col+"-+"
            head += "="*col+"=+"

        def show_contents(columns,items,withTrailing=True):
            if len(columns) != len(items):
                raise Exception("Dimension mismatch")

            toPrint = "| "

            for i,item in enumerate(items):
                item = str(item)

                if len(item) >= columns[i]:
                    raise Exception("col %s not large enough min = %s"%(i,len(item)+2))

                toPrint += item+" "*(columns[i]-len(item))+"| "

            print(toPrint[:-1])
            
            if withTrailing:
                print(row)

        print("\nAnnotation Summary\n---------------------")
        print(row)
        items = ['Taxon','genes (annots)','uniprot (annots)','coding-genes','combined','coverage']
        show_contents(columns,items,withTrailing=False)
        print(head)

        for taxon in self.taxaList:
            summary = self.get_annotation_summary(taxon,useIea=useIea)

            if float(summary["coverage"]) < minCoverage:
                continue

            items = [taxon,
                     "%s (%s)"%(summary["num_gene"],summary["num_annotated_genes"]),
                     "%s (%s)"%(summary["num_uniprot"],summary["num_annotated_uniprot"]),
                     summary["num_coding_genes"],
                     summary["num_annotated_combined"],
                     summary["coverage"]]
            show_contents(columns,items)

    def get_annotation_summary(self,taxonId,useIea=False):
        """
        return the annotation summary information for a taxa
        """
        
        acceptedCodes = get_evidence_codes(useIea=useIea)
        taxaQuery = self.session.query(Taxon).filter_by(ncbi_id=taxonId).first()
        annotations = self.session.query(GoAnnotation).filter(GoAnnotation.taxa_id==taxaQuery.id).\
                      filter(GoAnnotation.evidence_code.in_(acceptedCodes)).all()

        geneIds = [g.id for g in self.session.query(Gene).filter_by(taxa_id=taxaQuery.id).all()]
        uniprotQuery = self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()

        def remove_empty(lst):
            if None in lst:
                lst.remove(None)

        ## get how many genes code for proteins
        codingQuery = self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()
        codingGenes = list(set([u.gene_id for u in uniprotQuery]))
        remove_empty(codingGenes)
        
        ## get number of genes/proteins with at least one annotation
        annotatedGenes = list(set([a.gene_id for a in annotations]))
        annotatedProts = list(set([a.uniprot_id for a in annotations]))
        remove_empty(annotatedGenes)
        remove_empty(annotatedProts)

        ## remove genes covered that are not of the correct taxa (i.e. viral)
        apQuery = [self.session.query(Uniprot).filter_by(id=uid).first() for uid in annotatedProts]
        _genesFromUniprot = [self.session.query(Gene).filter_by(id=uq.gene_id).first() for uq in apQuery]
        
        if not _genesFromUniprot:
            geneIdsFromUniprot = []
        else:
            remove_empty(_genesFromUniprot)
            genesFromUniprot = []

            for gene in _genesFromUniprot:
                if gene == None:
                    continue
                geneTaxonId = self.session.query(Taxon).filter_by(id=gene.taxa_id).first()
                if int(geneTaxonId.ncbi_id) == int(taxonId):
                    genesFromUniprot.append(gene)
          
            geneIdsFromUniprot = [str(g.id) for g in genesFromUniprot]
        
        annotatedGenes = [str(a) for a in annotatedGenes]

        ## get total genes covered (uniprot + ncbi)
        genesCovered = list(set(annotatedGenes).union(set(geneIdsFromUniprot)))

        ## calculate coverage
        if len(genesCovered) == 0 or len(codingGenes) == 0:
            percentage = 0.0
        else:
            percentage = '%s'%(round((float(len(genesCovered))/float(len(geneIds))) * 100.0,4))

        results = {"num_gene":len(geneIds),
                   "num_uniprot":len(uniprotQuery),
                   "num_coding_genes":len(codingGenes),
                   "num_annotated_genes":len(annotatedGenes),
                   "num_annotated_uniprot":len(annotatedProts),
                   "num_annotated_combined":len(genesCovered),
                   "coverage":percentage}

        return results

    def read_file(self,taxaFilePath):
        """
        Use csv to read in the species 
        """
        
        fid = open(taxaFilePath,'rU')
        reader = csv.reader(fid,delimiter='\t')
        taxaList = []

        for linja in reader:
            if re.search("\D",linja[0]):
                continue
            taxaList.append(linja[0])

        fid.close()

        return taxaList

if __name__ == "__main__":
    print "Running..."

    tsum = TaxaSummary(["8364"])
    tsum.print_names_summary()
    tsum.print_annotation_summary()
