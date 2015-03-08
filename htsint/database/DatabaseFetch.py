#!/usr/bin/env python
"""
fetch all data files for htsint
At least 20GB of free space is recommended
"""
                                                          
import os,sys,subprocess,re,time,csv
from sys import platform as _platform
from htsint import __basedir__
from htsint import Configure

class DatabaseFetch(object):
    """
    Run each time we want files updated for database
    """

    def __init__(self,wget=os.path.join("/","usr","bin","wget"),gunzip=os.path.join("/","usr","bin","gunzip")):
        """  
        Constructor
        """
 
        config = Configure()
       
        ## check that we are in linux or osx
        if _platform == "linux" or _platform == "linux2":
            pass
        elif _platform == "darwin":
            pass
        elif _platform == "win32":
            raise Exception("DatabaseFetch currently does not work on windows platforms\n"+\
                            "you may still download the files one by one and populated the db")

        ## ensure they have set up config
        for key in ['data','dbname']:
            if config.log[key] == '':
                raise Exception("You must modify the config file before running DatabaseFetch.py")

        dataDir = config.log['data']

        if not os.path.isdir(dataDir):
            raise Exception("Specified htsint data directory does not exist %s"%dataDir)

        ## move into specified directory
        self.cwd = os.getcwd()
        os.chdir(dataDir)

        ## check for the necessary programs
        for wgetPath in [wget,os.path.join('/','usr','local','bin','wget')]:
            if os.path.exists(wgetPath) == True:
                break
            wgetPath = None

        if wgetPath == None:
            raise Exception("ERROR: cannot find wget -- either download the files (see documentation) or specify a path")
        self.wgetPath = wgetPath
        
        gunzipPaths = [gunzip,os.path.join("/","bin","gunzip")]
        for gunzipPath in gunzipPaths:
            if os.path.exists(gunzipPath) == True:
                break
            gunzipPath = None

        if gunzipPath == None:
            raise Exception("ERROR: cannot find gunzip -- either download the files (see documentation) or specify a path")
        self.gunzipPath = gunzipPath

    def _run_subprocess(self,cmd):
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stdin=subprocess.PIPE)
        while True:
            try:
                next_line = proc.stdout.readline()
                proc.wait()
                if next_line == '' and proc.poll() != None:
                    break
            except:
                proc.wait()
                break

    def fetch_file(self,fetchURL):
        cmd = "%s -N %s"%(self.wgetPath,fetchURL)
        self._run_subprocess(cmd)
        
    def unzip_file(self,fileName):
        print('unzipping...%s'%fileName)
        cmd = "%s -c %s > %s.db"%(self.gunzipPath,fileName,fileName[:-3])
        self._run_subprocess(cmd)

    def untar_file(self,fileName):
        print('untarring...%s'%fileName)
        cmd = "tar xzf %s"%fileName
        self._run_subprocess(cmd)

    def run(self):
        ## prepare a log file
        fid = open('fetchdb.log','w')
        writer = csv.writer(fid)

        def push_out(line):
            writer.writerow([line])
            print(line)

        push_out(sys.argv[0])
        push_out(time.asctime())
        push_out("fetching files...")

        ## fetch the files required for the database
        uniprotUrl = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/"
        filesToFetch = ["ftp://ftp.geneontology.org/pub/go/ontology/go.obo",
                        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
                        "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",
                        "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
                        "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz",
                        uniprotUrl + "idmapping/idmapping.dat.gz",
                        uniprotUrl + "idmapping/LICENSE",
                        uniprotUrl + "complete/uniprot_sprot.fasta.gz"
                    ]
                                 
        for fetchURL in filesToFetch:
            fileName = os.path.split(fetchURL)[-1]
            push_out("fetching %s..."%fileName)
            timeStart = time.time()
            self.fetch_file(fetchURL)
            fetchTime = time.time() - timeStart
            push_out("...%s"%fetchTime)

            ## unzip the gz files
            if not re.search("\.gz",fileName):
                continue

            if not os.path.exists(fileName[:-3]+".db") or fetchTime > 10:
                if re.search("\.tar\.gz",fileName):
                    self.untar_file(fileName)
                else:
                    self.unzip_file(fileName)

        push_out("complete.")
        fid.close()

        ## move back to original directory
        os.chdir(self.cwd)


if __name__ == "__main__":
    print "Running..."
    dbf = DatabaseFetch()
    dbf.run()
