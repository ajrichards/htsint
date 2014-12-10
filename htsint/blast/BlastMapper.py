#!/usr/bin/env python
"""
a class to handle results from BLAST searches

the parsed file that BlastMapper handles is a csv file with the following rows
    query,hit-identifier,hit-identifier-long,e-score,bit-score

    query is specific to the assembly
    HitId is a RefSeq

there is an example of this file in the unittest directory
It can be created with ParseBlast.py or ParallelParseBlast.py

"""

import os,sys,csv,re,getopt,time
import numpy as np
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene,Uniprot,Refseq
from htsint.database import get_idmapping_file, gene_mapper,taxa_mapper,uniprot_mapper
import gc

__author__ = "Adam Richards"

class BlastMapper(object):
    """
    Takes a parsed blast output and wraps it into a convenience class
    """

    def __init__(self):
        """
        Constructor
        """

        self.session,self.engine = db_connect()
        self.conn = self.engine.connect()

    def get_gene_id(self,geneDbId,gene2id):
        """
        returns the gene and species from uniprot
        """
        if str(geneDbId) == 'None':
            return '-','-'

        geneId,speciesId = '-','-'         
        for _speciesId, geneDict in gene2id.iteritems():
            if geneDict.has_key(geneDbId) and geneDict[geneDbId] != 'None':
                geneId = geneDict[geneDbId]
                speciesId = _speciesId
        
        return geneId,speciesId

    def create_summarized(self,parsedFilePath,summaryFilePath=None):
        """
        htsint uses output 5 (XML) and then parses it into a simple csv file
        """
        
        ## error checking
        if not os.path.exists(parsedFilePath):
            raise Exception("cannot find parsed file")

        if summaryFilePath == None:
            summaryFilePath = re.sub("\.csv","",parsedFilePath)+ "_summary.csv"

        ## input/output
        fidin = open(parsedFilePath,'rU')
        reader = csv.reader(fidin)
        header = reader.next()

        ## prepare out file         
        fidout = open(summaryFilePath,'w')
        writer = csv.writer(fidout)
        writer.writerow(["queryId","hitId","hitNcbiId","hitSpecies","hitSpeciesNcbiId","e-value"])

        ## read through the blast file and extract a unique list of ids 
        uniprotEntries = set([])
        for linja in reader:
            hitIdLong = linja[2]
            _hitId  = hitIdLong.split(" ")[1].split("|")
            if _hitId[-1] == '':
                hitId = _hitId[-2]
            else:
                hitId = _hitId[-1]
                
            uniprotEntries.update([hitId])
                 
        uniprotEntries = list(uniprotEntries)
        fidin.close()

        print("batch querying %s UniProt entries in the database... this may take some time"%(len(uniprotEntries)))
        #results = self.conn.execute(mytable.__table__.select(mytable.value.in_(values))
        #vailable_values = set(row.value for row in results)
        results = self.conn.execute(Uniprot.__table__.select(Uniprot.uniprot_entry.in_(uniprotEntries)))
        upEntry2Gene = dict([(str(row.uniprot_entry),str(row.gene_id)) for row  in results])
        upEntry2Taxa = dict([(str(row.uniprot_entry),str(row.taxa_id)) for row  in results])

        #s = select([Uniprot.id,Uniprot.uniprot_entry,Uniprot.gene_id,Uniprot.taxa_id]).where(Uniprot.uniprot_entry.in_(uniprotEntries))
        #_upQueries = self.conn.execute(s)
        #upQueries = _upQueries.fetchall()
        #upEntry2Gene = dict([(str(uquery['uniprot_entry']),str(uquery['gene_id'])) for uquery in upQueries])
        #upEntry2Taxa = dict([(str(uquery['uniprot_entry']),str(uquery['taxa_id'])) for uquery in upQueries])

        for key, item in upEntry2Gene.iteritems():
              print key, item
        
        ## query the taxa just to double check
        taxaList = list(set(upEntry2Taxa.values()))
        s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.name.in_(taxaList))
        _taxaQueries = self.conn.execute(s)
        taxaQueries = _taxaQueries.fetchall()
        taxaId2Ncbi = dict([(str(tquery['id']),str(tquery['ncbi_id'])) for tquery in taxaQueries])
        taxaName2Ncbi = dict([(str(tquery['name']),str(tquery['ncbi_id'])) for tquery in taxaQueries])
        for key, item in taxaName2Ncbi.iteritems():
            print key, item

        ## create a single dictionary of all gene information 
        gene2id = {}
        for tquery in taxaQueries:
            s = select([Gene.id,Gene.ncbi_id],Gene.taxa_id==tquery['id'])
            _geneQueries = self.conn.execute(s)
            taxaDict = dict([(str(r['id']),str(r['ncbi_id'])) for r in _geneQueries.fetchall()])
            print("there are  %s genes from %s (%s)"%(len(taxaDict.keys()),tquery['name'],tquery['ncbi_id']))
            gene2id[str(tquery['ncbi_id'])] = taxaDict

        ## read through the file again
        fidin = open(parsedFilePath,'rU')
        reader = csv.reader(fidin)
        header = reader.next()

        for linja in reader:
            query = linja[0]
            hitIdShort = linja[1]
            hitIdLong = linja[2]
            eScore = linja[3]
            bitScore = linja[4]
            queryId = query.split(" ")[0]            

            _hitId  = hitIdLong.split(" ")[1].split("|")
            if _hitId[-1] == '':
                hitId = _hitId[-2]
            else:
                hitId = _hitId[-1]
            
            hitNcbiId,hitSpeciesNcbiId = '-','-'
            hitSpecies = re.findall("OS=.+GN=",hitIdLong)[0][3:-4]
            
            geneId,hitSpeciesNcbiId = self.get_gene_id(hitId,gene2id)
            writer.writerow([queryId,hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,eScore])

        fidin.close()
        fidout.close()

    def load_summary_file(self,filePath,taxaList=None):
        """
        Because BLAST searches can result in large results lists there are several options
        
        filePath - path to summaryFile
        taxaList - if specified ignore any taxa not found in list otherwise take all taxa
        
        At this point the function returns only the best match for a each taxa

        """

        ## error checking
        if not os.path.exists(filePath):
            raise Exception("cannot find summary file")

        ## if taxa filter provided
        if taxaList != None:
            s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(taxaList))
            _taxaQueries = self.conn.execute(s)
            taxaQueries = _taxaQueries.fetchall()
            selectedTaxa = [str(tquery['id']) for tquery in taxaQueries]
            taxa2name = dict([(str(tquery['id']),str(tquery['ncbi_id'])) for tquery in taxaQueries])

        ## read the results
        results = {}
        fid = open(resultsFilePath,'rU')
        reader = csv.reader(fid)
        header = reader.next()

        uniqueQueries = set([])
        totalQueries = 0
        unfilteredQueries = 0

        for linja in reader:
            print linja
            if len(linja) == 4:
                queryId = linja[0]
                hitId = linja[1]
                hitNcbiId = linja[2]
                _evalue = float(linja[3])
            else:
                queryId = linja[0]
                queryNcbi = linja[1]
                hitId = linja[2]
                hitNcbiId = linja[3]
                _evalue = linja[4]

        if asGenes == True:
            queryId = re.sub("_i\d+","",queryId)




if __name__ == "__main__":
    print "Running..."
