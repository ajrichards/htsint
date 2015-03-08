#!/usr/bin/env python
"""
These are the tables for the database
"""

### make imports
import sys,os,re
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence
Base = declarative_base()

class Taxon(Base):
    '''
    class that handles the ncbi taxa
    '''

    __tablename__ = 'taxa'

    id = Column(Integer, Sequence('taxon_id_seq'),primary_key=True)
    ncbi_id = Column(Integer)
    name = Column(String)
    common_name_1 = Column(String)
    common_name_2 = Column(String)
    common_name_3 = Column(String)
    
    def __init__(self,ncbi_id,name='',common_name_1='',common_name_2='',common_name_3=''):
        self.ncbi_id = ncbi_id
        self.name = name
        self.common_name_1 = common_name_1
        self.common_name_2 = common_name_2
        self.common_name_3 = common_name_3
            
    def __repr__(self):
        return "<Taxon('%s','%s','%s','%s','%s')>"%(self.ncbi_id,
                                                    self.name,
                                                    self.common_name_1,
                                                    self.common_name_2,
                                                    self.common_name_3)

def taxa_mapper(session,ncbiIdList=None,myDict={},yieldPer=5):
    """
    a function that maps ncbi_ids to taxa.id
    if a dict is provided keys must be the string of the ncbi_id
    """

    if ncbiIdList:
        ncbiIdList = dict([(i,None)for i in list(set(ncbiIdList))])

    for t in session.query(Taxon).yield_per(yieldPer):
        if ncbiIdList and not ncbiIdList.has_key(str(t)):
            continue

        myDict[str(t.ncbi_id)] = t.id

    return myDict


class Gene(Base):
    '''
    class that handles ncbi genes
    '''

    __tablename__ = 'genes'

    id = Column(Integer, Sequence('gene_id_seq'),primary_key=True)
    ncbi_id = Column(String)
    description = Column(String,nullable=False)
    symbol = Column(String)
    synonyms = Column(String)
    taxa_id = Column(Integer,ForeignKey('taxa.id'),nullable=True)
    
    def __init__(self,ncbi_id,description,symbol,synonyms,taxa_id):
        self.ncbi_id = ncbi_id
        self.description = description
        self.symbol = symbol
        self.synonyms = synonyms
        self.taxa_id = taxa_id
               
    def __repr__(self):
        return "<Gene('%s','%s','%s','%s','%s')>"%(self.ncbi_id,
                                                   self.description,
                                                   self.symbol,
                                                   self.synonyms,
                                                   self.taxa_id)

def gene_mapper(session,ncbiIdList=None,myDict={}):
    """
    a function that maps ncbi_ids to gene.id
    if a dict is provided keys must be the string of the ncbi_id
    """

    if ncbiIdList:
        ncbiIdList = dict([(i,None)for i in list(set(ncbiIdList))])

    for g in session.query(Gene).yield_per(5):
        if ncbiIdList and not ncbiIdList.has_key(g):
            continue
        myDict[str(g.ncbi_id)] = g.id

    return myDict

class Uniprot(Base):
    '''
    class that handles Uniprot Accessions
    '''

    __tablename__ = 'uniprot'

    id = Column(Integer, Sequence('uniprot_id_seq'),primary_key=True)
    uniprot_ac = Column(String)
    uniprot_entry = Column(String)
    refseq = Column(String)
    taxa_id = Column(Integer,ForeignKey('taxa.id'),nullable=True)
    gene_id = Column(Integer,ForeignKey('genes.id'),nullable=True)

    def __init__(self,uniprot_ac,uniprot_entry,refseq,taxa_id,gene_id):
        self.uniprot_ac = uniprot_ac
        self.uniprot_entry = uniprot_entry
        self.refseq = refseq
        self.taxa_id = taxa_id
        self.gene_id = gene_id

    def __repr__(self):
        return "<Uniprot('%s','%s','%s','%s','%s')>"%(self.uniprot_ac,
                                                      self.uniprot_entry,
                                                      self.refseq,
                                                      self.taxa_id,
                                                      self.gene_id)

def uniprot_mapper(session,uniprotIdList=None,myDict={},gene=False,taxa=False):
    """
    a function that maps uniprot_entry to uniprot.id
    if a dict is provided keys must be the string of the ncbi_id

    session - current sqlalchemy session
    uniprotIdList - returns all entries in the table unless a list of uniProt_entries are provided
    myDict - provides a means to expand or update an existing dictionary
    gene - appends to the returned values the database gene id
    taxa - appends to the returned values the database taxa_id

    """

    if uniprotIdList:
        uniprotIdList = dict([(i,None) for i in list(set(uniprotIdList))])

    for u in session.query(Uniprot).yield_per(5):
        if uniprotIdList and not uniprotIdList.has_key(u.uniprot_entry):
            continue

        if gene == False and taxa == False:
            myDict[str(u.uniprot_entry)] = u.id
        elif gene == True and taxa == False:
            myDict[str(u.uniprot_entry)] = {'id':u.id,'gene_id':u.gene_id}
        elif gene == False and taxa == True:
            myDict[str(u.uniprot_entry)] = {'id':u.id,'taxa_id':u.taxa_id}
        elif gene == True and taxa == True:
            myDict[str(u.uniprot_entry)] = {'id':u.id,'gene_id':u.gene_id,'taxa_id':u.taxa_id}

    return myDict

class GoTerm(Base):
    '''
    class that handles gene ontology terms
    '''

    __tablename__ = 'go_terms'

    id = Column(Integer, Sequence('go_term_id_seq'),primary_key=True)
    go_id = Column(String)
    aspect = Column(String)
    name = Column(String)
    alternate_id = Column(String)
    description = Column(String)

    def __init__(self,go_id,aspect,name,alternate_id,description):
        self.go_id = go_id
        self.aspect = aspect
        self.name = name
        self.alternate_id = alternate_id
        self.description = description
    
    def __repr__(self):
        return "<GoTerm('%s','%s','%s','%s','%s')>"%(self.go_id,
                                                     self.aspect,
                                                     self.name,
                                                     self.alternate_id,
                                                     self.description)

def goterm_mapper(session,gotermIdList=None,myDict={}):
    """
    a function that maps uniprot_id to uniprot.id
    if a dict is provided keys must be the string of the ncbi_id
    """

    if gotermIdList:
        gotermIdList = dict([(i,None)for i in list(set(gotermIdList))])

    for g in session.query(GoTerm).yield_per(5):
        if gotermIdList and not gotermIdList.has_key(g):
            continue

        myDict[g.go_id] = g.id

    return myDict

class GoAnnotation(Base):
    '''
    class that handles gene ontology annotations
    annotations may be associated with a gene id or a uniprot id
    '''

    __tablename__ = 'go_annotations'

    id = Column(Integer, Sequence('go_annotation_id_seq'),primary_key=True)
    go_term_id = Column(Integer, ForeignKey('go_terms.id'))
    evidence_code = Column(String)
    pubmed_refs = Column(String)
    uniprot_id = Column(Integer, ForeignKey('uniprot.id'),nullable=True)
    gene_id = Column(Integer, ForeignKey('genes.id'),nullable=True)
    taxa_id = Column(Integer,ForeignKey('taxa.id'),nullable=False)

    def __init__(self,go_term_id,evidence_code,pubmed_refs,uniprot_id,gene_id,taxa_id):
        self.go_term_id = go_term_id
        self.evidence_code = evidence_code
        self.pubmed_refs = pubmed_refs
        self.uniprot_id = uniprot_id
        self.gene_id = gene_id
        self.taxa_id = taxa_id

    def __repr__(self):
        return "<GoAnnotation('%s','%s','%s','%s','%s','%s')>"%(self.go_term_id,
                                                                self.evidence_code,
                                                                self.pubmed_refs,
                                                                self.uniprot_id,
                                                                self.gene_id,
                                                                self.taxa_id)
