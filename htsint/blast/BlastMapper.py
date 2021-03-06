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

    def create_summarized(self,parsedFilePath,summaryFilePath=None,large=False,uniprot=False,species=None,
                          taxaList=[],hit2gene=None):
        """
        htsint uses output 5 (XML) and then parses it into a simple csv file
        large - use True if the parsed file as more than a few hundred hits
        uniprot - use True if the target database used is a Uniprot database
        species - if None htsint will try to find the species otherwise the scientific name is given
        taxaList - assign when BLAST is against a db other than uniprot
        hit2gene - a dictionary of keys that match the hits with values that match ncbi gene ids
        """

        ## error checking
        if not os.path.exists(parsedFilePath):
            raise Exception("cannot find parsed file")

        if not uniprot and len(taxaList) == 0:
            raise Exception("databases other than Uniprot must have a taxaList")    
        
        if summaryFilePath == None:
            summaryFilePath = re.sub("\.csv","",parsedFilePath)+ "_summary.csv"

        ## input/output
        timeStart = time.time()
        fidin = open(parsedFilePath,'r')
        reader = csv.reader(fidin)
        header = reader.__next__()

        ## prepare out file         
        fidout = open(summaryFilePath,'w')
        writer = csv.writer(fidout)
        writer.writerow(["queryId","hitId","hitNcbiId","hitSpecies","hitSpeciesNcbiId","e-value"])

        ## read through the blast file and extract a unique list of ids 
        hitEntries = set([])
        for linja in reader:
            hitIdLong = linja[2]
            _hitId  = hitIdLong.split(" ")[1].split("|")
            if _hitId[-1] == '':
                hitId = _hitId[-2]
            else:
                hitId = _hitId[-1]

            hitEntries.update([hitId])

        fidin.close()    
        hitEntries = list(hitEntries)

        if uniprot:
            print("batch querying %s UniProt entries in the database... this may take some time"%(len(hitEntries)))
            upEntry2Gene, upEntry2Taxa = {},{}
            if large == False:
                results = self.conn.execute(Uniprot.__table__.select(Uniprot.uniprot_entry.in_(hitEntries)))
                for row in results:
                    upEntry2Gene[str(row.uniprot_entry)] = str(row.gene_id)
                    upEntry2Taxa[str(row.uniprot_entry)] = str(row.taxa_id)
            else:
                ## using htsint's mapper
                uMapper = uniprot_mapper(self.session,uniprotIdList=hitEntries,gene=True,taxa=True)
                for key,item in uMapper.items():
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

        ## read through the file again
        fidin = open(parsedFilePath,'r')
        reader = csv.reader(fidin)
        header = reader.__next__()

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

            ## extract species id
            if species:
                hitSpecies = species
            elif re.search("OS=.+[A-Z]=",hitIdLong):
                hitSpecies = re.findall("OS=.+[A-Z]=",hitIdLong)[0][3:-4]
            elif re.search("\[.+\]",hitIdLong):
                hitSpecies = re.findall("\[.+\]",hitIdLong)[0][1:-1]
            else:
                print("WARNING: cannot find hitSpecies\n%s"%hitIdLong)
                hitSpecies = '-'
            
            ## clean species name    
            if re.findall("[A-Z]=",hitSpecies):
                hitSpecies = hitSpecies[:re.search("[A-Z]=",hitSpecies).start()-2]

            if uniprot and hitId in upEntry2Gene and str(upEntry2Gene[hitId]) != 'None':
                if upEntry2Gene[hitId] in gene2id and str(gene2id[upEntry2Gene[hitId]]) != 'None':
                    hitNcbiId = gene2id[upEntry2Gene[hitId]]

            if hit2gene and hitId in hit2gene:
                hitNcbiId = hit2gene[hitId]

            ## get taxa id associated with uniprot id
            if uniprot:
                hitSpeciesNcbiId = '-'
                if hitId in upEntry2Taxa and str(upEntry2Taxa[hitId]) != 'None':
                    if upEntry2Taxa[hitId] in taxaId2Ncbi and str(upEntry2Taxa[hitId]) != 'None':
                        hitSpeciesNcbiId = taxaId2Ncbi[upEntry2Taxa[hitId]]
                    else:
                        print('Were in! but taxaId2Ncbi does not have', upEntry2Taxa[hitId])
            elif len(taxaList) == 1:
                hitSpeciesNcbiId = taxaList[0]
                        
            writer.writerow([queryId,hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,eScore])

        print("blast summarize: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
        fidin.close()
        fidout.close()
        
        return summaryFilePath

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
        fid = open(filePath,'r')
        reader = csv.reader(fid)
        header = reader.__next__()

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
            if _evalue > evalue:
                evalueFilter += 1
                continue

            if taxaList and hitSpeciesNcbiId not in taxaList:
                taxaFilter += 1
                continue

            ## use the best evalue 
            if queryId not in results:
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

        fid.close()            
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
        for transcriptId, items in bmap.items():
            hitId,hitNcbiId,hitSpecies,hitSpeciesNcbiId,evalue = items
    
            if removeStrain:
                hitSpecies = re.sub("\s+$","",re.sub("\(.+\)","",hitSpecies))

            if hitSpecies not in taxaHits:
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
        for ln1,ln2 in longNames.items():
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

        print('total species:', len(taxaHits.keys()))
        print('check', np.array(sizes).sum())

        ## save results as csv file
        if not csvName:
            return

        fid = open(csvName,'w')
        writer = csv.writer(fid)
        writer.writerow(["SpeciesName","TotalHits","Percentage"])
        for i in range(taxaNames.size):
            writer.writerow([taxaNames[i],taxaCounts[i],percents[i]])

        fid.close()

    def print_summary(self,bmap):
        """
        print a bmap summary
        queries are usually transcripts
        hits depend on the database
        """

        hits,genes = set([]),set([])
        for key,item in bmap.iteritems():
            if type(item) == type([]):
                hitIds = [i[0] for i in item]
                geneIds = [i[1] for i in item]
            else:
                hitIds = [item[0]]
                geneIds = [item[1]]

            hits.update(hitIds)
            genes.update(geneIds)
        hits = list(hits)
        genes = list(genes)
        if '-' in genes:
            genes.remove("-")
        
        print("...queries: %s"%(len(bmap.keys())))
        print("...genes: %s"%(len(genes)))
        print("...hits: %s"%(len(hits)))
        return genes

    def get_gene_dict(self,bmap):
        """
        returns a gene to transcript dictionary
        """
         
        gene2transcript = {}
        for key,item in bmap.items():

            if type(item) == type([]):
                genes =  [i[1] for i in item]
            else:
                genes = [item[1]]

            for gene in genes:
                if gene == '-':
                    continue
                if gene not in gene2transcript:
                    gene2transcript[gene] = []
                gene2transcript[gene].append(key)

        return gene2transcript


    def get_hit_dict(self,bmap):
        """
        returns a hit to transcript dictionary
        """
         
        hit2transcript = {}
        for key,item in bmap.items():

            if type(item) == type([]):
                hits =  [i[0] for i in item]
            else:
                hits = [item[0]]

            for hit in hits:
                if hit == '-':
                    continue
                if hit not in hit2transcript:
                    hit2transcript[hit] = []
                hit2transcript[hit].append(key)

        return hit2transcript



if __name__ == "__main__":
    print("Running...")
