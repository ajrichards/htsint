.. main file for lpedit documentation

Installation
======================

hts-integrate (``htsint``) is intended to run under GNU/Linux or another UNIX-based environment.  The software is tested with Ubuntu 14.04 and OSX although it should run under any operating system as it is simply a Python package with an attached PostgreSQL (or another) database.  Because Debian based distros are common in the bioinformatics world, some installation details are Debian specific, though they may be modified for other systems.  

Running local BLAST through hts-integrate is optional.  For BLAST functionality of the package install the Suite of `BLAST+ tools <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_. 

   .. code-block:: bash
     
      ~$ sudo apt-get install ncbi-blast+  


Install hst-integrate
______________________

``htsint`` is at its core a Python package so it may be installed (or upgraded) like most other packages in the `Python package index <Python Package Index>`_.


   .. code-block:: bash

      ~$ pip install htsint

   .. code-block:: bash

      ~$ pip install htsint --upgrade


.. note:: hts-integrate can also be installed from source.

Setting up your database
_______________________________

   1. Obtain the source for ``htsint``

      .. code-block:: bash

         ~$ git clone https://github.com/ajrichards/hts-integrate.git
         

   2. Create a `PostgreSQL <http://www.postgresql.org>`_ database.

      If PostgreSQL is not installed.  See the `postgres docs <http://www.postgresql.org/docs>`_ for further detail.

      .. code-block:: bash
   
         ~$ sudo apt-get install postgresql
       
      Create a database (the ``postgres=#`` indicates the postgres prompt)

      .. code-block:: bash

         ~$ sudo su - postgres
         ~$ psql -U postgres
         postgres=# CREATE USER htsint WITH ENCRYPTED PASSWORD '*********';
         postgres=# CREATE DATABASE htsintdb WITH OWNER hstint;
         postgres=# \q

   3. Configure ``htsint`` to see the database

      Create a config file.

      .. code-block:: bash

         ~$ cd hts-integrate/htsint
         ~$ cp configure.py.dist configure.py

      Then edit the file so that it looks something like this.

      .. code-block:: python

          CONFIG = {'data':'/usr/local/share/htsint',
          'dbname':"htstintdb",
          'dbuser':"htsint",
          'dbpass':"*********",
          'dbhost':"localhost",
          'dbport':"5432"
          }
 
      If you leave the ``dbpass`` field blank ``htsint`` will prompt you each time it needs access to the database.
      Also, be sure that you change the ``data`` directory to somewhere ``htsint`` can read and write.


Installing from source
___________________________________

Installation from source requires that the prerequsites are present before proceeding.

First install the prerequsite Python packages.

      .. code-block:: bash

         ~$ sudo apt-get install python-numpy python-matplotlib python-networkx python-biopython
         ~$ sudo apt-get install python-scipy python-pymc

Then install the database (see previous section for details on setup).

      .. code-block:: bash

         ~$ sudo apt-get install python-psycopg2 python-sqlalchemy 
        
Install ``htsint``

      .. code-block:: bash

         ~$ cd /hts-integrate
         ~$ sudo python setup.py install

Unittests
^^^^^^^^^^^^^^

First you will need to download the data and populate the database.  Instructions for this are on the :doc:`database <database>` page.  Once everything is setup you can run the unittests to ensure everything is working properly.

   .. code-block:: bash 
      
         ~$ cd /hts-integrate
         ~$ python runTests.py

Links
__________
 
   * `How to run BLAST on a local computer <http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local>`_
   * `Pip documentation <https://pip.readthedocs.org/en/latest/>`_
