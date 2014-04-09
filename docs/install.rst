.. main file for lpedit documentation

Installation
======================

``htsint`` is a Python package.  The NGS analysis pipelines available through ``htsint`` make extensive use of the Python ecosystem so there are a fair number of prerequsites.  ``htsint`` is developed and tested under Linux and OSX operating systems so the installation is specific to these systems, however, all of the prerequisites are operating system independent and these instructions could be modified for Windows and other operating systems fairly easily.


Setting up your system
--------------------------

These are the instructions that are tested under Ubuntu 13.10.  They may be modified for other distros or operating systems by following the links and installing as specified for your OS.

   1. Obtain the source for ``htsint``

      .. code-block: bash

         ~$ git co https://github.com/ajrichards/hts-integrate.git
         

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

   4. Install the prerequsites and package itself


      The core packages.

      .. code-block:: bash

         ~$ sudo apt-get install python-numpy python-matplotlib python-networkx python-biopython

      The database packages.

      .. code-block:: bash

         ~$ sudo apt-get install python-psycopg2 python-sqlalchemy ncbi-blast+

      The statistics packages.

      .. code-block:: bash
 
         ~$ sudo apt-get install python-scipy python-pymc

      Install ``htsint``

      .. code-block:: bash

         ~$ cd /hts-integrate
         ~$ sudo python setup.py install


Next you will need to download the data and populate the database.  Instructions for this are on the :doc:`database <database>` page.  Once everything is setup you can run the unittests to ensure everything is working properly.

   .. code-block:: bash 
      
         ~$ cd /hts-integrate
         ~$ python runTests.py
