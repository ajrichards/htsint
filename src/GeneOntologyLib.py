#!/usr/bin/env python
"""
library of functions for use with the GeneOntology class
"""

import os,re
from htsint import __basedir__

def get_ontology_file():
    """
    check for presence of ontology file
    raise exception when not found
    return the file path
    """

    ontologyFile = os.path.join(__basedir__,'database','go.obo')
    if os.path.exists(ontologyFile) == False:
        raise Exception("Could not find 'go.obo' -- did you run FetchGo.py?")

    return ontologyFile


def read_ontology_file():
    """
    read the ontology file to find term-term edges
    """

    ontologyFile = get_ontology_file()

    fid = open(ontologyFile,'r')
    termCount = 0
    for linja in fid.readlines():
        linja = linja[:-1]
        print linja
        if re.search("^\[Term\]",linja):
            termCount += 1
        


        if termCount == 2:
            break
