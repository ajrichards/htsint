#!/usr/bin/python 
import unittest,sys,os,getopt

from htsint import __basedir__
sys.path.append(__basedir__)

import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('agg')

from unittests import *
unittest.main()
