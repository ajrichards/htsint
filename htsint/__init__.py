import sys,os

## basic files
#sys.path.append(os.path.realpath(__file__))
#print(os.path.realpath(__file__))
#print(os.path.pardir(os.path.realpath(__file__)))
#print(os.path.split(os.path.abspath(__file__))[:-1])

from .version import __version__
from .basedir import __basedir__

## database functions and classes
from .Configure import Configure
from .RunSubprocess import run_subprocess, RunSubprocess
from .GeneOntology import GeneOntology
from .TermDistances import TermDistances
from .GeneDistances import GeneDistances
from .AssembleDistances import AssembleDistances
from .TaxaSummary import TaxaSummary 
from .GeneSetCollection import GeneSetCollection
from .GeneSet import GeneSet
