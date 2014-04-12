#!/usr/bin/env python

import sys
import couchdb
from htsint.database import Couch
from htsint import __basedir__
sys.path.append(__basedir__)
try:
    from configure import CONFIG
except:
    CONFIG = {'dbhost':'localhost','dbport':'5984'}

#class HtsDb(Couch):
#    """
#    A generic class
#    """


class HtsDb(object):

    def __init__(self):
        """ 
        initialize class to handle several couchdbs 
        """

        self.dbnames = ['htsint-uniprot','htsint-gene','htsint-go','htsint-taxa']
        self.db = {}
        self.server = couchdb.Server("http://%s:%s/"%(CONFIG['dbhost'],CONFIG['dbport']))
        
    def connect(self):
        """
        connect to the htsint specific databases
        """

        for dbname in self.dbnames:
            self.db[dbname] = self.server[dbname]

    def save_doc(self,dbname,doc):
        """
        save a document into a database
        """

        self.db[dbname].save(doc)

# If your CouchDB server is running elsewhere, set it up like this:
#couch = couchdb.Server('http://example.com:5984/')

# select database
#db = couch['mydb']

#create a document and insert it into the db:
#doc = {'foo': 'bar'}
#db.save(doc)
