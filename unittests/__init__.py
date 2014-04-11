import unittest,getopt,sys,os
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('agg')

## parse inputs                                                                                                                      
try:
    optlist, args = getopt.getopt(sys.argv[1:],'v')
except getopt.GetoptError:
    print getopt.GetoptError
    print sys.argv[0] + "-v"
    print "... the verbose flag (-v) may be used"
    sys.exit()

try:
    from configure import CONFIG
except:
    CONFIG = None

if CONFIG == None:
    raise Exception("You must create a configure.py before running the unittests")

VERBOSE = False
RUNALL = False

for o, a in optlist:
    if o == '-v':
        VERBOSE = True

## Database tests
from DatabaseTest import *
DatabaseTestSuite = unittest.TestLoader().loadTestsFromTestCase(DatabaseTest)
DatabaseSuite = unittest.TestSuite([DatabaseTestSuite])

## GeneOntology tests
from GeneOntologyTest import *
GeneOntologyTestSuite = unittest.TestLoader().loadTestsFromTestCase(GeneOntologyTest)
GeneOntologySuite = unittest.TestSuite([GeneOntologyTestSuite])

## Blast Tests
if CONFIG['blast'] == True:
    from BlastTest import *
    BlastTestSuite = unittest.TestLoader().loadTestsFromTestCase(BlastTest)
    BlastSuite = unittest.TestSuite([BlastTestSuite])
