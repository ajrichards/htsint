#!/usr/bin/env python
"""
library of functions for use with the GeneOntology class
"""

import os,sys,re
from htsint import __basedir__

sys.path.append(__basedir__)
try:
    from configure import CONFIG
except:
    CONFIG = None

def get_ontology_file():
    """
    check for presence of ontology file
    raise exception when not found
    return the file path
    """

    if CONFIG == None:
        raise Exception("You must create a configure.py before GeneOntology")

    dataDir = CONFIG['data']

    ontologyFile = os.path.join(dataDir,'go.obo')
    if os.path.exists(ontologyFile) == False:
        raise Exception("Could not find 'go.obo' -- did you run FetchGo.py?")

    return ontologyFile

def read_ontology_file():
    """
    read the ontology file to find term-term edges
    store all the relationships in a dictionary
    """

    ontologyFile = get_ontology_file()
    fid = open(ontologyFile,'r')
    termCount = 0
    goId = None
    goDict = {"cellular_component":{},
              "molecular_function":{},
              "biological_process":{}}

    def add_term(goNamespace,source,sink):
        #print source,sink
        if len(sink) != 1 or source == sink[0]:
            return
        if not re.search("GO\:",source) or not re.search("GO\:",sink[0]):
            raise Exception("Invalid go id in ontology file: %s, %s"%(source,sink[0]))

        if goDict[goNamespace].has_key(source) == False:
            goDict[goNamespace][source] = set([])
        goDict[goNamespace][source].update(sink)

    for linja in fid.readlines():
        linja = linja[:-1]

        ## find go id
        if re.search("^id\:",linja):
            goId = re.sub("^id\:|\s+","",linja)
            goName,goNamespace = None,None
            is_a,part_of = None,None
            isObsolete = False
            termCount += 1
            continue 

        ## find a few go id attributes
        if re.search("^name\:",linja):
            goName = re.sub("^name\:\s+","",linja)
        if re.search("^namespace\:",linja):
            goNamespace = re.sub("^namespace\:\s+","",linja)
        if re.search("^def\:",linja):
            goDef = re.sub("^def\:\s+","",linja)
            if re.search("OBSOLETE\.",goDef):
                isObsolete = True

        ## ignore obolete terms
        if goId == None or isObsolete == True:
            continue

        ## add term relationships (is_a)
        if re.search("^is\_a\:",linja):
            is_a = re.sub("^is\_a\:\s+","",linja)
            is_a_sink = re.findall("GO\:\d+",is_a)
            add_term(goNamespace,goId,is_a_sink)

        ## other term relationships
        #if re.search("^relationship\:",linja):
        #    part_of_sink = re.findall("GO\:\d+",linja)
        #    add_term(goNamespace,goId,part_of_sink)
            
    return goDict
