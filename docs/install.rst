.. main file for lpedit documentation

Installation
======================

To install ``htsint`` there are two steps.

   1. install the python package
   2. setup the database

``htsint`` is tested with current versions of Ubuntu and OSX although it should run under any operating system as it is simply a Python package with an attached PostgreSQL database.  Because Debian based distros are common in the bioinformatics world, some installation details are Debian specific.  

To ensure that you can run BLAST locally install the Suite of `BLAST+ tools <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_. For example,

   .. code-block:: bash
     
      ~$ sudo apt-get install ncbi-blast+

As of the current release the following Python packages are used by ``htsint``.

   * `NumPy <http://www.numpy.org>`_
   * `SciPy <http://www.scipy.org>`_
   * `matplotlib <http://matplotlib.org>`_
   * `NetworkX <https://networkx.github.io>`_
   * `Biopython <http://biopython.org/wiki/Main_Page>`_
   * `Psycopg2 <http://initd.org/psycopg>`_
   * `SQLAlchemy <http://www.sqlalchemy.org>`_

Install hstint
-------------------

``htsint`` is at its core a Python package so it may be installed (or upgraded) like most other packages in the `Python package index <Python Package Index>`_.

   .. code-block:: bash

      ~$ pip install htsint

   .. code-block:: bash

      ~$ pip install htsint --upgrade


.. note:: ``htsint`` can also be installed from source.

Setting up your database
---------------------------

Follow the :doc:`database setup instructions <database>`. Ensure that you edit the ``taxa`` field of the config file to fit your needs.

Installing from source
-----------------------------

Installation from source requires that the prerequsites are present before proceeding.

First install the prerequsite Python packages.  For example,

   .. code-block:: bash

      ~$ sudo apt-get install python-numpy python-matplotlib python-networkx python-biopython
      ~$ sudo apt-get install python-scipy python-psycopg2 python-sqlalchemy 

Then install the database (see previous section for details on setup).
        
Install ``htsint``

   .. code-block:: bash

      ~$ cd /htsint
      ~$ sudo python setup.py install

Unittests
^^^^^^^^^^^^^^

First you will need to download the data and populate the database.  Instructions for this are on the :doc:`database <database>` page.  Once everything is setup you can run the unittests to ensure everything is working properly.

   .. code-block:: bash 
      
      ~$ cd /htsint
      ~$ python runTests.py

Links
----------
 
   * `How to run BLAST on a local computer <http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local>`_
   * `Pip documentation <https://pip.readthedocs.org/en/latest/>`_
