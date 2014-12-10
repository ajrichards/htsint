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

    def get_gene_id(self,hitId,uniprot2id,gene2id):
        """
        returns the gene and species from uniprot
        """
        geneId,speciesId = '-','-'         
        for _speciesId, uniprotDict in uniprot2id.iteritems():
            if uniprotDict.has_key(hitId) and uniprotDict[hitId] != 'None':
                if gene2id.has_key(uniprotDict[hitId]):
                    geneId = gene2id[uniprotDict[hitId]]
                speciesId = _speciesId
        if speciesId == '-':
            print hitId, [(sid,uniprotDict.__contains__(hitId)) for sid,uniprotDict in uniprot2id.iteritems()]

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
        writer.writerow(["queryId","hitId","hitNcbi","hitSpecies","e-value"])

        ## read through the blast file and extract the taxa
        taxaList = set([])
        for linja in reader:
            hitIdLong = linja[2]
            hitSpecies = re.findall("OS=.+GN=",hitIdLong)[0][3:-4]
            taxaList.update([hitSpecies])

        taxaList = list(taxaList)
        print(taxaList)
        fidin.close()

        # query all the relevant taxa
        s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.name.in_(taxaList))
        _taxaQueries = self.conn.execute(s)
        taxaQueries = _taxaQueries.fetchall()
        taxa2name = dict([(str(tquery['id']),str(tquery['ncbi_id'])) for tquery in taxaQueries])

        ## create a single dictionary of all gene information 
        gene2id = {}
        for tquery in taxaQueries:
            s = select([Gene.id,Gene.ncbi_id],Gene.taxa_id==tquery['id'])
            _geneQueries = self.conn.execute(s)
            taxaDict = dict([(str(r['id']),str(r['ncbi_id'])) for r in _geneQueries.fetchall()])
            print("there are  %s genes from %s (%s)"%(len(taxaDict.keys()),tquery['name'],tquery['ncbi_id']))
            gene2id.update(taxaDict)

        ## create a dictionary of uniprot names (per species)
        uniprot2id = {}

        for tquery in taxaQueries:
            s = select([Uniprot.uniprot_entry,Uniprot.gene_id],Uniprot.taxa_id==tquery['id'])
            _uniprotQueries = self.conn.execute(s)
            taxaDict = dict([(str(r['uniprot_entry']),str(r['gene_id'])) for r in _uniprotQueries.fetchall()])
            print("there are  %s uniprot from %s (%s)"%(len(taxaDict.keys()),tquery['name'],tquery['ncbi_id']))
            uniprot2id[str(tquery['ncbi_id'])] = taxaDict
                
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
            hitNcbiId,queryNcbiId = '-','-'
            _hitId  = hitIdLong.split(" ")[1].split("|")

            ## id fetch fix
            if _hitId[-1] == '':
                hitId = _hitId[-2]
            else:
                hitId = _hitId[-1]

            queryId = query.split(" ")[0]
            
            geneId,speciesId = self.get_gene_id(hitId,uniprot2id,gene2id)

            #if speciesId == '-':
            #    taxa2name

            writer.writerow([queryId,hitId,geneId,speciesId,eScore])

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
