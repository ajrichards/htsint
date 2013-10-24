'''
Created Oct 24th, 2013
'''

import os,sys
from distutils.core import setup
from distutils.extension import Extension
from numpy import get_include

sys.path.append("src")
from version import __version__

__author__ = "AJ Richards"

from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
try:
    from py2exe.build_exe import py2exe
except:
    pass

sys.path.append("lpEdit")
from version import __version__

DESCRIPTION = "Cross-platform editor to facilitate literate programming"
LONG_DESCRIPTION = """                                                                                                                                                                                            

A package meant to construct Bayesian models that aid in the discovering of
biomarkers via the integration high-throughput data sources with additional 
covariates.

"""

REQUIRES = ['numpy', 'matplotlib','pymc']
DISTNAME = 'lpEdit'
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

ISRELEASED = True
VERSION = __version__
FULLVERSION = VERSION
if not ISRELEASED:
    FULLVERSION += '.beta'

def configuration(parent_package='', top_path=None):
    config = Configuration(None, parent_package, top_path,
                           version=FULLVERSION)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('src')    
    config.add_data_dir(os.path.join('src','examples'))

    return config

if __name__ == '__main__':
    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=AUTHOR_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          long_description=LONG_DESCRIPTION,
          classifiers=CLASSIFIERS,
          platforms='any',
          configuration=configuration)
