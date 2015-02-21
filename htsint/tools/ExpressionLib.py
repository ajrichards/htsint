#!/usr/bin/env python

"""
tools for expression and count based tasks

"""

import os,sys,csv,gc,re
import numpy as np

def read_RSEM_counts_files(geneFilePath,isoformFilePath):
    """
    read the RSEM counts files into a matrix
    """

    if not os.path.exists(geneFilePath):
        raise Exception("Cannot find gene file\n%s"%(geneFilePath))
    if not os.path.exists(isoformFilePath):
        raise Exception("Cannot find isoform file\n%s"%(isoformFilePath))

    ## load the gene counts
    fid1 = open(geneFilePath,'rU')
    reader1 = csv.reader(fid1,delimiter="\t")
    header1 = reader1.next()    
    results1 = {}
    check = 0
    gc.disable()
    
    for linja in reader1:
        check += 1
        results1[linja[0]] = {'transcript':linja[1],'length':float(linja[2]),'eff_length':float(linja[3]),\
                              'exp_count':int(round(float(linja[4]))),'TPM':float(linja[5]),'FPKM':float(linja[6])}
    fid1.close()
    if check != len(results1.keys()):
        raise Exception("Rows in gene count file are not first columns unique")

    ## load the isoform results
    fid2 = open(isoformFilePath,'rU')
    reader2 = csv.reader(fid2,delimiter="\t")
    header2 = reader2.next()    
    results2 = {}
    check = 0

    for linja in reader2:
        check += 1
        results2[linja[0]] = {'gene':linja[1],'length':float(linja[2]),'eff_length':float(linja[3]),\
                              'exp_count':float(linja[4]),'TPM':float(linja[5]),'FPKM':float(linja[6])}
    fid1.close()
    if check != len(results2.keys()):
        raise Exception("Rows in gene count file are not first columns unique")

    fid2.close()
    gc.enable()

    return results1, results2

def read_matrix(matFilePath,delimiter=",",mtype='float'):
    """
    assumes that row one are the samples and col one are the transcripts
    matrix can only be of mtype 'int' or 'float'

    """

    print 'reading', matFilePath

    if mtype not in ['int','float']:
        raise Exception("mtype must be 'int' or 'float'")
    if not os.path.exists(matFilePath):
        raise Exception("Cannot find matFilePath\n%s"%matFilePath)

    fid = open(matFilePath,'rb')
    reader = csv.reader(fid,delimiter=delimiter)
    header = reader.next()

    ## get the gene and sample ids
    transcriptIds = []
    sampleIds = np.array(header[1:])
    gc.disable()
    for linja in reader:
        transcriptIds.append(linja[0])
    gc.enable()
    transcriptIds = np.array(transcriptIds)
    fid.close()

    ## fill in the matrix
    mat = np.zeros((transcriptIds.shape[0],sampleIds.shape[0]),dtype=mtype)
    fid = open(matFilePath,'rb')
    reader = csv.reader(fid,delimiter=delimiter)
    header = reader.next()
    row = 0 
    for linja in reader:
        if mtype == 'int':
            mat[row,:] = [int(float(i)) for i in linja[1:]]
        else:
            mat[row,:] = [float(i) for i in linja[1:]]
        row +=1
    fid.close()

    return transcriptIds,sampleIds,mat

def read_de_results(filePath,delimiter=",",tool="edgeR"):
    """
    read the differential expression output from DESeq or edgeR

    """

    print 'reading', filePath

    if not os.path.exists(filePath):
        raise Exception("Cannot find matFilePath\n%s"%filePath)

    if tool not in ["edgeR","DESeq"]:
        raise Exception("invalid tool specified use 'edgeR' or 'DESeq'")

    fid = open(filePath,'rb')
    reader = csv.reader(fid,delimiter=delimiter)
    
    ## get columnIds
    header = reader.next()
    columnIds = np.array(header[1:])

    ## get the gene and sample ids
    transcriptIds = []

    gc.disable()
    for linja in reader:
        transcriptIds.append(linja[0])
    gc.enable()
    transcriptIds = np.array(transcriptIds)
    fid.close()

    ## fill in the matrix
    mat = np.zeros((transcriptIds.shape[0],columnIds.shape[0]))
    fid = open(filePath,'rb')
    reader = csv.reader(fid,delimiter=delimiter)
    header = reader.next()
    row = 0 
    for linja in reader:
        _row = [re.sub("NA","NaN",i) for i in linja[1:]]             
        mat[row,:] = [float(i) for i in _row]
        row +=1
    fid.close()

    return transcriptIds,columnIds,mat
    [(x, y) for x in [1,2,3] for y in [3,1,4] if x != y]#

def create_count_matrix(results,label,sampleList):
    """
    this function is untested
    """

    ## use first sample to get rows
    mat = np.zeros((len(results[0].keys()),len(sampleList)))
    keys = sorted(np.array(results[0].keys()))

    for j,sample in enumerate(sampleList):
        for i,key in enumerate(keys):
            mat[i,j] = results[j][key]['exp_count']

    ## write to file 
    fid = open("%s-counts.csv"%label,'w')
    writer = csv.writer(fid)
    if re.search("gene",label):
        writer.writerow(["gene"]+sampleList)
    else:
        writer.writerow(["isoform"]+sampleList)

    for r in range(mat.shape[0]):
        row = [keys[r]] + [int(i) for i in mat[r,:].tolist()]
        writer.writerow(row)

    fid.close()
