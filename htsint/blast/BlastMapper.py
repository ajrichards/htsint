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
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene,Uniprot,uniprot_mapper

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
        self.hits = None

    def create_summarized(self,parsedFilePath,summaryFilePath=None,large=False):
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

        timeStart = time.time()
        upEntry2Gene, upEntry2Taxa = {},{}
        if large == False:
            results = self.conn.execute(Uniprot.__table__.select(Uniprot.uniprot_entry.in_(uniprotEntries)))
            for row in results:
                upEntry2Gene[str(row.uniprot_entry)] = str(row.gene_id)
                upEntry2Taxa[str(row.uniprot_entry)] = str(row.taxa_id)
        else:
            ## using htsint's mapper
            uMapper = uniprot_mapper(self.session,uniprotIdList=uniprotEntries,gene=True,taxa=True)
            for key,item in uMapper.iteritems():
                upEntry2Gene[str(key)] = str(item['gene_id'])
                upEntry2Taxa[str(key)] = str(item['taxa_id'])
       
        ## query the taxa just to double check
        taxaList = list(set(upEntry2Taxa.values()))
        while 'None' in taxaList: taxaList.remove('None')
        s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.id.in_([int(tid) for tid in taxaList]))
        _taxaQueries = self.conn.execute(s)
        taxaQueries = _taxaQueries.fetchall()
        taxaId2Ncbi = dict([(str(tquery['id']),str(tquery['ncbi_id'])) for tquery in taxaQueries])

        ## create a single dictionary of all gene information 
        gene2id = {}
        for taxaDbId in taxaList:
            s = select([Gene.id,Gene.ncbi_id],Gene.taxa_id==taxaDbId)
            _geneQueries = self.conn.execute(s)
            taxaDict = dict([(str(r['id']),str(r['ncbi_id'])) for r in _geneQueries.fetchall()])
            #print("there are  %s genes from %s (%s)"%(len(taxaDict.keys()),tquery['name'],tquery['ncbi_id']))
            gene2id.update(taxaDict)
        print("uniprot batch query: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))

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
            hitSpecies = re.findall("OS=.+[A-Z]=",hitIdLong)[0][3:-4]
            if re.findall("[A-Z]=",hitSpecies):
                hitSpecies = hitSpecies[:re.search("[A-Z]=",hitSpecies).start()-2]

            ## map the uniprot id to gene
            if upEntry2Gene.has_key(hitId) and str(upEntry2Gene[hitId]) != 'None':
                if gene2id.has_key(upEntry2Gene[hitId]) and str(gene2id[upEntry2Gene[hitId]]) != 'None':
                    hitNcbiId = gene2id[upEntry2Gene[hitId]]
       
            ## get taxa id associated with uniprot id
            hitSpeciesNcbiId = '-'
            if upEntry2Taxa.has_key(hitId) and str(upEntry2Taxa[hitId]) != 'None':
                if taxaId2Ncbi.has_key(upEntry2Taxa[hitId]) and str(taxaId2Ncbi.has_key(upEntry2Taxa[hitId])) != 'None':
                    hitSpeciesNcbiId = taxaId2Ncbi[upEntry2Taxa[hitId]]
                else:
                    print 'Were in! but taxaId2Ncbi does not have', upEntry2Taxa[hitId]
            else:
                print 'upEntry2Taxa does not have', hitId
            writer.writerow([queryId,hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,eScore])

        fidin.close()
        fidout.close()

    def load_summary(self,filePath,taxaList=None,trinityGene=False,evalue=0.0001,best=True):
        """
        Because BLAST searches can result in large results lists there are several options
        
        filePath - path to summaryFile
        taxaList - if specified ignore any taxa not found in list otherwise take all taxa
        
        At this point the function returns only the best match for a each taxa

        A summary results file has the following header:
        queryId,hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,e-value

        trinityGene (default = False) - function assumes that each queryId is associated with a unique identifier
                                           if set to True using Trinity nameing scheme results become 'gene' centric
        
        best (default = True) - if specified as True then only the query results with the lowest e-value is returned
                                otherwise all results are returned

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
        fid = open(filePath,'rU')
        reader = csv.reader(fid)
        header = reader.next()

        uniqueQueries = set([])
        totalQueries = 0
        totalHits = 0
        evalueFilter = 0
        taxaFilter = 0
        for linja in reader:
            if len(linja) == 6:
                queryId = linja[0]
                hitId = linja[1]
                hitNcbiId = linja[2]
                hitSpecies = linja[3]
                hitSpeciesNcbiId = linja[4]
                _evalue = float(linja[5])
            else:
                raise Exception("Invalid number of columns in line in summary blast file \n%s"%linja)

            if trinityGene == True:
                queryId = re.sub("_i\d+","",queryId)

            totalHits += 1

            # filtering
            if '-' in linja:
                continue
            if _evalue > evalue:
                evalueFilter += 1
                continue

            if taxaList and hitSpeciesNcbiId not in taxaList:
                taxaFilter += 1
                continue

            ## use the best evalue 
            if not results.has_key(queryId):
                if best:
                    results[queryId] = (hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,_evalue)
                else:
                    results[queryId] = [(hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,_evalue)]

            else:
                ## if not in append mode and we have smaller evalue
                if best and _evalue < results[queryId][4]:
                    results[queryId] = (hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,_evalue)
                ## if append mode and we have smaller evalue
                elif (hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,_evalue) in results[queryId]:
                    continue
                elif not best and _evalue < results[queryId][0][4]:
                    results[queryId] = [(hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,_evalue)] + results[queryId]
                elif not best and _evalue >= results[queryId][0][4]:
                    results[queryId].append((hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,_evalue))

        print("queries filtered due to evalue > %s: %s"%(evalue,evalueFilter))
        print("queries filtered due to taxa: %s"%(taxaFilter))
        print("total hits: %s"%totalHits)
        self.hits = totalHits
        
        return results


    def make_taxa_pie_chart_and_table(self,bmap,removeStrain=False,threshold=2,figName=None,csvName=None,
                                      ax=None,shadow=True,startangle=90,explode=True,fontSize=11):
        """
        summarize the taxa using a pie chart and reStructuredText table
        
        all species that make up less than threshold percent will be grouped as 'other'

        """

        taxaHits = {}
        for transcriptId, items in bmap.iteritems():
            hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,evalue = items
    
            if removeStrain:
                hitSpecies = re.sub("\s+$","",re.sub("\(.+\)","",hitSpecies))

            if not taxaHits.has_key(hitSpecies):
                taxaHits[hitSpecies] = 0
            taxaHits[hitSpecies] += 1

        taxaNames = np.array(taxaHits.keys())
        taxaCounts = np.array(taxaHits.values()).astype(float)
        percents = (taxaCounts / taxaCounts.sum()) * 100.0

        ## reorder as ranked
        rankedInds = np.argsort(taxaCounts)[::-1]
        taxaNames = taxaNames[rankedInds]
        taxaCounts = taxaCounts[rankedInds]
        percents = percents[rankedInds]

        ## create pie plot
        includedInds = np.where(percents >= threshold)[0]
        otherPercent = percents[np.where(percents < threshold)[0]].sum()
        
        if otherPercent > 0.0:
            labels = taxaNames[includedInds].tolist() + ['other']
            sizes = percents[includedInds].tolist() + [otherPercent]
        else:
            labels = taxaNames[includedInds].tolist()
            sizes = percents[includedInds].tolist()
       
        ## explode the biggest chunk
        if explode:
            explodeInd = sizes.index(max(sizes))
            explode = np.zeros(len(sizes))
            explode[explodeInd] = 0.1
        else:
            explode = None

        if ax == None:
            fig = plt.figure(figsize=(9,7))
            ax = fig.add_subplot(111)
     
        colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral','lightgray',
                  'khaki','honeydew','tomato','cornflowerblue','darkseagreen','darkorchid'] * 100

        ## shorten long names
        longNames = {"Schizosaccharomyces pombe":"S. pombe"}
        for ln1,ln2 in longNames.iteritems():
            if ln1 in labels:
                lnind = labels.index(ln1)
                labels[lnind] = ln2
            
        patches, texts, autotexts =ax.pie(sizes, labels=labels,explode=explode,colors=colors,
                                          autopct='%1.1f%%', shadow=shadow, startangle=startangle)        
        ax.axis('equal')
        
        ## adjust font sizes
        proptease = fm.FontProperties()
        proptease.set_size(fontSize)
        plt.setp(autotexts, fontproperties=proptease)
        plt.setp(texts, fontproperties=proptease)

        if figName:
            plt.savefig(figName)

        print 'total species:', len(taxaHits.keys())
        print 'check', np.array(sizes).sum()

        ## save results as csv file
        if not csvName:
            return

        fid = open(csvName,'w')
        writer = csv.writer(fid)
        writer.writerow(["SpeciesName","TotalHits","Percentage"])
        for i in range(taxaNames.size):
            writer.writerow([taxaNames[i],taxaCounts[i],percents[i]])

        fid.close()

    def get_dicts(self,bmap):
        """
        returns a gene to transcript and transcript to gene dictionary
        """
         
        transcript2gene = {}
        gene2transcript = {}
        for key,_item in bmap.iteritems():
            for item in _item:
                if item[1] != '-':
                    transcript2gene[key] = item[1]
                    if not gene2transcript.has_key(item[1]):
                        gene2transcript[item[1]] = {}
                    if not gene2transcript[item[1]].has_key(item[3]):
                        gene2transcript[item[1]][item[3]] = set([])
                    gene2transcript[item[1]][item[3]].update([key])
        return transcript2gene,gene2transcript  


if __name__ == "__main__":
    print "Running..."
