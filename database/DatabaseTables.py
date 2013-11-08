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
    chromosome = Column(String)
    map_location = Column(String)
    synonyms = Column(String)
    taxa_id = Column(Integer,ForeignKey('taxa.id'))
    
    def __init__(self,ncbi_id,description,symbol,chromosome,map_location,
                 synonyms,taxa_id):
        self.ncbi_id = ncbi_id
        self.description = description
        self.symbol = symbol
        self.chromosome = chromosome
        self.map_location = map_location
        self.synonyms = synonyms
        self.taxa_id = taxa_id
               
    def __repr__(self):
        return "Gene('%s','%s','%s','%s','%s','%s')>"%(self.ncbi_id,
                                                       self.description,
                                                       self.symbol,
                                                       self.chromosome,
                                                       self.map_location,
                                                       self.synonyms,
                                                       self.taxa_id)
class Accession(Base):
    '''
    class that handles ncbi genes
    '''

    __tablename__ = 'accessions'

    id = Column(Integer, Sequence('accession_id_seq'),primary_key=True)
    status = Column(String)
    nucleo_gi = Column(String)
    protein_gi = Column(String)
    genomic_nucleo_gi = Column(String)
    genomic_start = Column(String)
    genomic_stop = Column(String)
    orientation = Column(String)
    assembly = Column(String)
    gene_id = Column(Integer,ForeignKey('genes.id'))
    
    def __init__(self,status,nucleo_gi,protein_gi,genomic_nucleo_gi,genomic_start,
                 genomic_stop,orientation,assembly,gene_id):
        self.status = status
        self.nucleo_gi = nucleo_gi
        self.protein_gi = protein_gi
        self.genomic_nucleo_gi = genomic_nucleo_gi
        self.genomic_start = genomic_start
        self.genomic_stop = genomic_stop
        self.orientation = orientation
        self.assembly = assembly
        self.gene_id = gene_id
               
    def __repr__(self):
        return "Accession('%s','%s','%s','%s','%s','%s','%s','%s','%s')>"%(self.status,
                                                                           self.nucleo_gi,
                                                                           self.protein_gi,
                                                                           self.genomic_nucleo_gi,
                                                                           self.genomic_start,
                                                                           self.genomic_stop,
                                                                           self.orientation,
                                                                           self.assembly,
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
    gene_id = Column(Integer, ForeignKey('genes.id'))

    def __init__(self,go_term_id,evidence_code,pubmed_refs,gene_id):
        self.go_term_id = go_term_id
        self.evidence_code = evidence_code
        self.pubmed_refs = pubmed_refs
        self.gene_id = gene_id
        
    def __repr__(self):
        return "<GoTermInstance('%s','%s','%s','%s')>"%(self.go_term_id,
                                                        self.evidence_code,
                                                        self.pubmed_refs,
                                                        self.gene_id)
"""
class Pub(Base):
    '''
    class that handles publications
    '''

    __tablename__ = 'pubs'

    id = Column(Integer,primary_key=True)
    pmc_id = Column(String)
    title = Column(String)
    
    def __init__(self,pmc_id,title):
        self.pmc_id = pmc_id
        self.title = title
            
    def __repr__(self):
        return "<Pub('%s','%s')>"%(self.pmc_id,
                                   self.title)
"""
"""
class PubInstance(Base):
    '''
    class that handles publication instances
    '''

    __tablename__ = 'pub_instances'

    id = Column(String,primary_key=True)
    pub_id = Column(String,ForeignKey('pubs.id'))
    gene_id = Column(String,ForeignKey('genes.id'))
    
    def __init__(self,pmc_id,gene_id):
        self.pmc_id = pmc_id
        self.gene_id = gene_id
            
    def __repr__(self):
        return "<PubInstance('%s','%s')>"%(self.pmc_id,
                                           self.genes_id)
"""

"""
class SemDistance(Base):
    '''
    class that handles publication instances
    '''

    __tablename__ = 'pub_instances'

    id = Column(String,primary_key=True)
    gene_id_1 = Column(String,ForeignKey('genes.id'))
    gene_id_2 = Column(String,ForeignKey('genes.id'))
    method =  Column(String)

    def __init__(self,gene_id_1,gene_id_2,method):
        self.gene_id_1 = gene_id_1
        self.gene_id_2 = gene_id_2
        self.method = method
            
    def __repr__(self):
        return "<SemDistancePub('%s','%s','%s')>"%(self.gene_id_1,
                                                   self.gene_id_2,
                                                   self.method)
"""
