#!/usr/bin/python

'''
htsint
config file
'''

__author__ = "A Richards"

import os,csv,re,ast
from version import __version__


defaultCONFIG = {'data':'/usr/local/share/htsint',
                 'dbname':"",
                 'dbuser':"",
                 'dbpass':"",
                 'dbhost':"localhost",
                 'dbport':"5432",
                 'taxa': ['3702','4932','5476','7227','7955','8355','8364','9031','9606',\
                          '10090','10566','10116','28377']
}


class Configure(object):
    """
    class to handle the database config
    """

    def __init__(self):
        """
        constructor
        """

        self.logFileDir = os.path.join(os.path.expanduser('~'),".hts-integrate")
        if os.path.exists(self.logFileDir) == True and os.path.isdir(self.logFileDir) == False:
            os.remove(self.logFileDir)
        
        if os.path.isdir(self.logFileDir) == False:
            os.mkdir(self.logFileDir)

        self.logFilePath = os.path.join(self.logFileDir,"dbconfig.log")

        
        if os.path.exists(self.logFilePath) == False:
            fid = open(self.logFilePath,'w')
            writer = csv.writer(fid)
            
            for key,item in defaultCONFIG.iteritems():
                if item == None:
                    item = 'None'
                elif type(item) != type('i am a string'):
                    item = str(item)

                writer.writerow([key,item])
            fid.close()

        self.log = self.read_project_log(self.logFilePath)
        
    ## effectivly the only action necessary to save a project in its current state
    def save(self):
        writer = csv.writer(open(self.logFilePath,'w'))

        for key,item in self.log.iteritems():
            if item == None:
                item = 'None'
            elif type(item) != type('i am a string'):
                item = str(item)

            writer.writerow([key,item])
            
    ## reads the log file assciated with the current project and returns a dict
    def read_project_log(self,logPathName):
        
        if os.path.isfile(logPathName) == False:
            print "ERROR: invalid model logfile specified",logPathName
            return None
        else:
            logFileDict = {}
            reader = csv.reader(open(logPathName,'r'))
            for linja in reader:

                if re.search("\[|\{|None",str(linja[1])):
                    try:
                        linja[1] = ast.literal_eval(str(linja[1]))
                    except:
                        print 'ERROR: Logger -- string literal conversion failed', linja[1]

                logFileDict[linja[0]] = linja[1]

            return logFileDict
