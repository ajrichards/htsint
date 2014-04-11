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
    taxa_id = Column(Integer,ForeignKey('taxa.id'))
    
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

class Uniprot(Base):
    '''
    class that handles Uniprot Accessions
    '''

    __tablename__ = 'uniprot'

    id = Column(Integer, Sequence('uniprot_id_seq'),primary_key=True)
    uniprot_id = Column(String)
    uniprot_entry = Column(String)
    refseq = Column(String)
    uniprot_taxa_id = Column(String)
    gene_id = Column(Integer,ForeignKey('genes.id'))

    def __init__(self,uniprot_id,uniprot_entry,refseq,taxa_id,gene_id):
        self.uniprot_id = uniprot_id
        self.uniprot_entry = uniprot_entry
        self.refseq = refseq
        self.uniprot_taxa_id = taxa_id
        self.gene_id = gene_id

    def __repr__(self):
        return "<Uniprot('%s','%s','%s','%s','%s')>"%(self.uniprot_id,
                                                      self.uniprot_entry,
                                                      self.refseq,
                                                      self.uniprot_taxa_id,
                                                      self.gene_id)

class GoTerm(Base):
    '''
    class that handles gene ontology terms
    '''

    __tablename__ = 'go_terms'

    id = Column(Integer, Sequence('go_term_id_seq'),primary_key=True)
    go_id = Column(String)
    aspect = Column(String)
    description = Column(String)
    
    def __init__(self,go_id,aspect,description):
        self.go_id = go_id
        self.aspect = aspect
        self.description = description
    
    def __repr__(self):
        return "<GoTerm('%s','%s','%s')>"%(self.go_id,
                                           self.aspect,
                                           self.description)

class GoAnnotation(Base):
    '''
    class that handles gene ontology annotations
    '''

    __tablename__ = 'go_annotations'

    id = Column(Integer, Sequence('go_annotation_id_seq'),primary_key=True)
    go_term_id = Column(Integer, ForeignKey('go_terms.id'))
    evidence_code = Column(String)
    pubmed_refs = Column(String)
    uniprot_id = Column(Integer, ForeignKey('uniprot.id'))
    taxa_id = Column(Integer,ForeignKey('taxa.id'))

    def __init__(self,go_term_id,evidence_code,pubmed_refs,uniprot_id,taxa_id):
        self.go_term_id = go_term_id
        self.evidence_code = evidence_code
        self.pubmed_refs = pubmed_refs
        self.uniprot_id = uniprot_id
        self.taxa_id = taxa_id

    def __repr__(self):
        return "<GoAnnotation('%s','%s','%s','%s')>"%(self.go_term_id,
                                                      self.evidence_code,
                                                      self.pubmed_refs,
                                                      self.uniprot_id,
                                                      self.taxa_id)
