#!/usr/bin/python

import time,csv,re,sys
import numpy as np
from htsint.database import get_annotation_file, get_gene2go_file


annotationFile = get_annotation_file()
annotationFid = open(annotationFile,'rU')
annotsCount = 0
annotatedIds = {}
annots1,annots2 = 0,0

for record in annotationFid:
    record = record[:-1].split("\t")
    if record[0][0] == "!":
        continue
    if record[0] != 'UniProtKB':
        continue

    taxon = re.sub("taxon:","",record[12])
    if taxon == "" or re.search("\|",taxon):
        continue

    annots1 += 1

gene2goFile = get_gene2go_file()
gene2goFid = open(gene2goFile,'rU')
header = gene2goFid.next()

for record in gene2goFid:
    annots2 += 1

print("__________________")
print("Uniprot annotations: %s"%annots1)
print("Gene2go annotations: %s"%annots2)
print("Total Annotations: %s"%(annots1 + annots2))
