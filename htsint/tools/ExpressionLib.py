#!/usr/bin/env python

"""
tools for expression and count based tasks

"""

import os,sys,csv,gc


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
