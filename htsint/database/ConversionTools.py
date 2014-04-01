#!/usr/bin/python

import os,re
import numpy as np



def convert_gene_ids(geneList,target):
    """
    takes a list of gene ids (int) and returns the target field
    Normally, the database can be used, however this script is 
    used for example if a species has not yet been entered in the db.

    """

    if target not in ['taxid','symbol']:
        raise Exception("Invalid target")

    geneInfoFile = os.path.join(os.path.split(os.path.abspath(__file__))[0],"gene_info.db")
    if os.path.exists(geneInfoFile) == False:
        raise Exception("ERROR: cannot find gene info file")

    geneInfoFID = open(geneInfoFile,'rU')
    header = geneInfoFID.next()

    if type(geneList) != type(np.array([])):
        geneList = np.array(geneList)

    toReturn = np.array([None for i in range(len(geneList))])
    found = 0

    for record in geneInfoFID:
        record = record.rstrip("\n")
        record = record.split("\t")
    
        if re.search("^\#",record[0]) or len(record) != 15:
            continue

        geneID = int(record[1])

        #if geneID not in geneList:
        #    continue
        
        indx = np.where(geneList == geneID)[0]
        
        if len(indx) == 0:
            continue
        
        found += 1
        indx = indx[0]

        if target == 'taxid':
            taxID = int(record[0])
            toReturn[indx] = taxID
        if target == 'symbol':
            symbol = record[2]
            toReturn[indx] = symbol

        #LocusTag = record[3]
        #Synonyms = record[4]
        #dbXrefs = record[5]
        #chromosome = record[6]
        #map_location = record[7]
        #description = record[8]
        #type_of_gene = record[9]
        #Symbol_from_nomenclature_authority = record[10]
        #Full_name_from_nomenclature_authority = record[11]
        #Nomenclature_status = record[12]
        #otherDesignations = record[13]
        #Modification_date = record[14]

    toReturn = toReturn.tolist()
    print "convert_gene_ids: %s/%s genes found."%(found,len(toReturn))

    return toReturn


if __name__ == "__main__":


    geneList = [101738635, 692505, 101738735, 101744422, 101267788, 101745190, 101738329, 101739693, 101260649,100808326]
    
    taxList = convert_gene_ids(geneList,'taxid')
    print taxList
