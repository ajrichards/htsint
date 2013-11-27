import sys,os

## basic files
from version import __version__
from basedir import __basedir__

## database functions and classes
sys.path.append(os.path.join(__basedir__,'database'))
from DatabaseTables import Base,Taxon,Gene,Accession,GoTerm,GoAnnotation
from DatabaseTools import db_connect
from ParseNcbiXml import ParseNcbiXml
from DatabaseAppend import DatabaseAppend
from PrimeDatabase import PrimeDatabase

