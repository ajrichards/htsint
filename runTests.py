#!/usr/bin/python 
import unittest,sys,os,getopt

import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('agg')

from unittests import *
unittest.main()
