#!/usr/bin/env python

import os,sys,re
from distutils.sysconfig import get_python_lib
from shutil import copytree
from numpy.distutils.core import setup

try:
    from py2exe.build_exe import py2exe
except:
    pass

sys.path.append("htsint")
from version import __version__

__author__ = "AJ Richards"

DESCRIPTION = "API and database for high-throughput sequencing analyses."
LONG_DESCRIPTION = """
A package meant to construct Bayesian models that aid in the discovering of
biomarkers via the integration high-throughput data sources with additional 
covariates.

"""

def get_files(dirPath):
    notIncluded = ["\.pyc"]
    allFiles = []
    for fileName in os.listdir(dirPath):
        include = True
        for pat in notIncluded:
            if re.search(pat,fileName):
                include = False
        filePath = os.path.join(dirPath,fileName)
        if include == True and os.path.isfile(filePath):
            allFiles.append(os.path.realpath(filePath))
            
    return allFiles

REQUIRES = ['numpy', 'matplotlib','pymc']
DISTNAME = 'htsint'
LICENSE = 'BSD 3-Clause'
AUTHOR = "Adam J Richards"
AUTHOR_EMAIL = "adamricha@gmail.com"
URL = 'https://github.com/ajrichards/hts-integrate'
CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'License :: OSI Approved :: BSD 3-Clause'
]

DATADIRS = ["database","blast","stats","tools"]
FILES = {'': [os.path.join("htsint","*.py")],
         'database': [os.path.join("htsint","database","*.py")],
         'blast': [os.path.join("htsint","blast","*.py *.pl")],
         'stats': [os.path.join("htsint","stats","*.py")],
         'tools':[os.path.join("htsint","tools","*.py")]}

ISRELEASED = True
VERSION = __version__
FULLVERSION = VERSION
if not ISRELEASED:
    FULLVERSION += '.beta'


if __name__ == '__main__':
    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=AUTHOR_EMAIL,
          description=DESCRIPTION,
          version=FULLVERSION,
          license=LICENSE,
          url=URL,
          package_dir={'': 'htsint'},
          long_description=LONG_DESCRIPTION,
          classifiers=CLASSIFIERS,
          package_data= FILES,
          platforms='any')
