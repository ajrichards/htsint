#!/usr/bin/env python
"""
A class that prints and generates summary information about a species list
with respect to the information contained in the database
"""

__author__ = "Adam Richards"

import csv,re,sys
from htsint.database import db_connect,Taxon,Gene,Uniprot,GoAnnotation,fetch_annotations

class TaxaSummary(object):
    """
    Summarize the information in a taxa list
    """

    def __init__(self,taxaFilePath,delimiter="\t",upass=''):
        """
        The input file is a tab delimited file 
        The delimiter can also be specified
        The first column is the ncbi taxa id
        All additional columns are optional and user specified 
        """

        self.session,self.engine = db_connect(upass=upass)
        self.taxaList = list(set(self.read_file(taxaFilePath)))
        print("There are %s unique taxa in the list"%(len(self.taxaList)))

    def print_names_summary(self):
        """
        print to the screen in REstructuredText
        """

        ## get the taxa ids associated wtih the taxa list
        col1,col2,col3 = 10,38,37
        row = "+"+"-"*col1+"-+"+"-"*col2+"-+"+"-"*col3+"-+"
        head = "+"+"="*col1+"=+"+"="*col2+"=+"+"="*col3+"=+"

        def show_contents(item1,item2,item3):
            item1 = str(item1)
            item2 = str(item2)
            item3 = str(item3)

            if len(item1) >= col1:
                raise Exception("col1 not large enough min = %s"%(len(item1)+3))
            if len(item2) >= col2:
                raise Exception("col2 not large enough min = %s"%(len(item2)+3))
            if len(item3) >= col3:
                raise Exception("col3 not large enough min = %s"%(len(item3)+3))

            print("| "+item1+" "*(col1-len(item1))+"| "+item2+" "*(col2-len(item2))+"| "+item3+" "*(col3-len(item3))+"|")
            print(row)

        show_contents('Taxon ID','Scientific Name','Common Name')
        print(head)

        for taxonId in self.taxaList:            
            query = self.session.query(Taxon).filter_by(ncbi_id=taxonId).all()
            if len(query) > 1:
                raise Exception("Taxon id '%s' return multiple values from database"%taxonId)
            if query == None:
                raise Exception("Invalid taxon '%s' provided in list"%taxonId)
            
            show_contents(query[0].ncbi_id, query[0].name, query[0].common_name_1)

    def print_annotation_summary(self):
        """
        print a summary of the annotation information
        """
        
        columns = [10,10,10,14,17,19,10,10]
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

        print(row)
        items = ['Taxon','genes (annots)','uniprot (annots)','coding-genes','combined','coverage']
        show_contents(columns,items,withTrailing=False)
        print(head)

        for taxon in self.taxaList:
            summary = self.get_annotation_summary(taxon)
            items = [taxon,
                     "%s (%s)"%(summary["num_gene"],summary["num_annotated_genes"]),
                     "%s (%s)"%(summary["num_uniprot"],summary["num_annotated_uniprot"]),
                     summary["num_coding_genes"],
                     summary["num_annotated_combined"],
                     summary["coverage"]]
            show_contents(columns,items)

    def get_annotation_summary(self,taxonId):
        """
        return the annotation summary information for a taxa
        """

        taxaQuery = self.session.query(Taxon).filter_by(ncbi_id=taxonId).first()
        annotations = self.session.query(GoAnnotation).filter_by(taxa_id=taxaQuery.id).all()
        geneIds = [g.id for g in self.session.query(Gene).filter_by(taxa_id=taxaQuery.id).all()]
        uniprotQuery = self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()

        def remove_empty(lst):
            if None in lst:
                lst.remove(None)

        ## get how many genes code for proteins
        codingQuery = self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()
        print('There are %s gene ids'%len(geneIds))
        print('There are %s uniprot ids'%len(uniprotQuery))
        codingGenes = list(set([u.gene_id for u in uniprotQuery]))
        remove_empty(codingGenes)
        #print('There are %s uniprot ids that code for proteins'%(len(codingGenes)))
        
        ## get number of genes/proteins with at least one annotation
        annotatedGenes = list(set([a.gene_id for a in annotations]))
        annotatedProts = list(set([a.uniprot_id for a in annotations]))
        remove_empty(annotatedGenes)
        remove_empty(annotatedProts)

        apQuery = self.session.query(Uniprot).filter(Uniprot.id.in_(annotatedProts)).all()
        genesFromUniprot = list(set([uq.gene_id for uq in apQuery]))
        remove_empty(genesFromUniprot)
        genesCovered = list(set(annotatedGenes).union(set(genesFromUniprot)))

        #print("There are %s genes with at least one annotation "%(len(annotatedGenes)))
        #print("There are %s proteins with at least one annotation"%(len(annotatedProts)))
        #print("There are %s genes with at least one annotation (ncbi + uniprot info)"%(len(genesCovered)))

        ## calculate coverage
        if len(genesCovered) == 0 or len(codingGenes) == 0:
            percentage = 0.0
        percentage = '%s'%(round((float(len(genesCovered))/float(len(codingGenes))) * 100.0,4))
        #print('Coverage: %s/%s'%(len(genesCovered),len(codingGenes)) + ' ( %s'%percentage + " % )")

        results = {"num_gene":len(geneIds),
                   "num_uniprot":len(uniprotQuery),
                   "num_coding_genes":len(codingGenes),
                   "num_annotated_genes":len(annotatedGenes),
                   "num_annotated_uniprot":len(annotatedProts),
                   "num_annotated_combined":len(genesCovered),
                   "coverage":percentage}

        return results

    def foo(self):
        print(row)
        sys.exit()

        ## print taxa name information
        print("------------------")
        print("taxa names")
        print("------------------")
        for tr in taxaRows:
            print("id='%s',name='%s',common_name='%s'"%(tr.ncbi_id,tr.name,tr.common_name_1))
    
        ## print taxa specific infomation
        print("------------------")
        print("taxa specific info")
        print("------------------")
        annots = {}
        
        for taxaRow in taxaRows:
            queryGenes = self.session.query(Gene).filter_by(taxa_id=taxaRow.id).all()
            geneList = [g.ncbi_id for g in queryGenes]
            print("id='%s',num_genes='%s'"%(taxaRow.ncbi_id,len(queryGenes)))
            
            #uniprotQueries = session.query(Uniprot).filter(Uniprot.uniprot_id.in_(identifiers)).all()
            #queryUniprot = self.session.query(Uniprot).filter(Uniprot.gene_id.in_([g.id for g in queryGenes])).all()
            #print("id='%s',num_uniprot_ids_associated_with_gene_ids='%s'"%(taxaRow.ncbi_id,len(queryUniprot)))
            #print [qu.uniprot_id for qu in queryUniprot]

            #a = self.session.query(GoAnnotation).filter_by(taxa_id=taxaRow.id).all()
            #print("blah %s"%len(a))

            #sys.exit()

            ## annotations associated with geneid or uniprot id tied to gene id
            geneUniAnnots = fetch_annotations(geneList,self.session,idType='ncbi',aspect='biological_process')
            print("id='%s',num_geneuni_annots='%s'"%(taxaRow.ncbi_id,len(geneUniAnnots)))
            for key,item in geneUniAnnots.iteritems():
                print key,item



            #self.session.query(GoAnnotation).filter_by(taxa_id=taxaId).all()

            ## annotations associated with uniprot ids that have taxa
            #annots[taxaRow.ncbi_id] = {'bytaxa':self.session.query(GoAnnotation).filter_by(taxa_id=taxaId).all()}



        sys.exit()
        annots = {}

        for t,taxaId in enumerate(taxaIds):
            print t,taxaId
            annots[self.taxaList[t]] = {'bytaxa':self.session.query(GoAnnotation).filter_by(taxa_id=taxaId).all()}

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
