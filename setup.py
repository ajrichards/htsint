#!/usr/bin/env python

import os,sys,re
from distutils.sysconfig import get_python_lib
from shutil import copytree
from distutils.core import setup
#from numpy.distutils.core import setup

try:
    from py2exe.build_exe import py2exe
except:
    pass

sys.path.append("htsint")
from version import __version__

__author__ = "AJ Richards"

DESCRIPTION = "API and database for high-throughput sequencing analyses."
LONG_DESCRIPTION = """
A package that creates gene sets using unsupervised clustering. Gene sets are
generated based on Gene Ontology information from one or more species.

"""

REQUIRES = ['numpy', 'matplotlib','networkx','sqlalchemy','biopython']
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
    'License :: OSI Approved :: BSD License'
]


FILES = {'htsint': [os.path.join('database','*.py'),
                    os.path.join('stats','*.py'),
                    os.path.join('blast','*.py'),
                    os.path.join('tools','*.py')]}

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
          packages=['htsint'],
          long_description=LONG_DESCRIPTION,
          classifiers=CLASSIFIERS,
          package_data=FILES,
          platforms='any')
