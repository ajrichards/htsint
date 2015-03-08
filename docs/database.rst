.. main file for lpedit documentation

Database
======================

The schema for the database.

.. figure:: dbschema.png
   :scale: 99%
   :align: center
   :alt: database schema
   :figclass: align-center


The main uses of the database are:

   * :doc:`name mapping <name-mapping>`
   * :doc:`annotation fetching <annotation-fetching>`


Setting up your database
-----------------------------

There are 4 step that you need to carry out in order to have a working database.

1. Install `PostgreSQL <http://www.postgresql.org>`_ and create an empty database
2. Modify the `hts-integrate` config file
3. Fetch all necessary data files
4. Populate the datbase

(1) PostgreSQL setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install `PostgreSQL <http://www.postgresql.org>`_ on you system.  For example on Ubuntu 

   .. code-block:: bash

      ~$ sudo apt-get install postgresql

Then you need to create a database.  Here we do this with a database named `htsintegrate`, a user named `dbuser` and a password `somepassword`.  You should change these parameters to fit your system.

   .. code-block:: bash

      ~$ sudo su - postgres
      ~$ psql -U postgres
      CREATE USER dbuser WITH ENCRYPTED PASSWORD 'somepassword';
      CREATE DATABASE htsintegrate WITH OWNER dbuser; 		   
      \q

Visit the `PostgreSQL documentation <http://www.postgresql.org/docs>`_ to learn more.

(2) Modify the config file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set the database name, user, password and port

>>> from htsint import Configure
>>> config = Configure()
>>> config.log['dbuser'] = 'dbuser'
>>> config.log['dbname'] = 'htsintegrate'
>>> config.log['dbpass'] = 'somepassword'
>>> config.save()

The other variables than can be modified are:

>>> print(config.log.keys())
['dbport', 'taxa', 'dbhost', 'dbpass', 'dbuser', 'blast', 'data', 'dbname']

Check that all taxa that you would like annotation information for are present.  You can print the default taxa with

>>> print(config.log['taxa'])
['10566', '9606', '10090', '8364', '10116', '9913', '8355', '9031', '9601', '7955', '9823', '9986', '9615',\
 '9541', '9598', '7227', '10029', '100141', '9796', '3702']

If you would like to add *Solanum lycopersicum* then use typical list syntax.  You need to save any changes you make.

>>> config.log['taxa'].append('4081')
>>> config.save()

It may be useful to specify the directory where all the downloaded data is stored (use full path)

The default is

>>> print(config.log['data'])
'/usr/local/share/htsint'

Alternatively, it is possible to edit this file directly.  To locate the directory where it is stored you may type the following.

>>> import os
>>> os.path.join(os.path.expanduser('~'),".hts-integrate")
'/home/adam/.hts-integrate'

The dbport (default '5432'), dbhost (default 'localhost') may also be configured.

.. note:: hts-integrate will only populate annotation information for taxa in the *taxa* variable so make sure all species are present **before** database population.


(3) Fetch the necessary data files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The rest of the process is automated assuming you have created your PostgreSQL database and modified you config file.

>>> from htsint.database import DatabaseFetch
>>> fetch = DatabaseFetch()
>>> fetch.run()

This class only currently works under Linux/OSX operating systems.  For other systems the following files could be downloaded by hand and placed in the 'data' directory.

   * `go.obo <ftp://ftp.geneontology.org/pub/go/ontology/go.obo>`_
   * `taxdump.tar.gz <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_
   * `gene_info.gz <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz>`_
   * `gene2go.gz <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz>`_
   * `gene_association.goa_uniprot.gz <ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz>`_
   * `idmapping.dat.gz <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping.dat.gz>`_
   * `uniprot_sprot.fasta.gz <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz>`_

A logfile is produced and stored in your data directory.

(4) Populate the database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the database can be populated with the following class.

>>> from htsint.database DatabaseCreate
>>> db = DatabaseCreate()
>>> db.run()

A logfile is produced and stored in your data directory.

Additional Notes
-----------------

Database updating
^^^^^^^^^^^^^^^^^^^^^^

Because of the challenges that can arise through naming conflicts when updating NCBI and UniProt data it is recommended that you run the fetch and create steps again, which will create a clean updated version.  The fetch step will check if a current file is the newest and only download a new one if necessary.

Database portability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may want to copy your database to another computer instead of waiting for the db to populate.  This can be done as follows.

To create a file that may be transferred to another computer

   .. code-block:: bash

      ~$ pg_dump -h localhost -U dbuser dbname > htsint.sql

To add the database to another server

   .. code-block:: bash

      ~$ sudo su - postgres
      ~$ psql -U postgres
      CREATE USER dbuser WITH ENCRYPTED PASSWORD 'somepassword';
      CREATE DATABASE newdbname WITH OWNER dbuser; 		   
      \q
      ~$ psql newdbname < htsint.sql
