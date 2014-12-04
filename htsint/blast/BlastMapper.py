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
        
        pass


    def load_parsed(self,parsedFilePath):
        """
        htsint uses output 5 (XML) and then parses it into a simple csv file
        """
        
        ## error checking
        if not os.path.exists(parsedFilePath):
            raise Exception("cannot find parsed file")

        fidin = open(parsedFilePath,'rU')
        reader = csv.reader(fidin)
        header = reader.next()

        ## read through the blast file and extract the taxa
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

            print 'query', query
            print '_hitId', _hitId
            print 'hitId', hitId
            print 'hitIdLong', hitIdLong
            print 'hitIdshort', hitIdShort

            queryId = query.split(" ")[0]
            hitSpecies = re.findall("OS=.+GN=",hitIdLong)[0][:-4]
            print hitSpecies
            sys.exit()

            #hitNcbiId = self.refseq2ncbi(hitId)

            sys.exit()

        fidin.close()

if __name__ == "__main__":
    print "Running..."
