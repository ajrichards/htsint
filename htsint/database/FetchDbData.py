#!/usr/bin/env python
"""
fetch all data files for htsint
At least 20GB of free space is recommended
"""
                                                          
import os,sys,subprocess,re,time,csv
from htsint import __basedir__
sys.path.append(__basedir__)

try:
    from configure import CONFIG
except:
    CONFIG = None

if CONFIG == None:
    raise Exception("You must create a configure.py before running FetchDbData.py")

dataDir = CONFIG['data']

if not os.path.isdir(dataDir):
    raise Exception("Specified htsint data directory does not exist %s"%dataDir)

## move into specified directory
cwd = os.getcwd()
os.chdir(dataDir)

## check for the necessary programs
for wgetPath in [os.path.join("/","usr","bin","wget"), os.path.join('/','usr','local','bin','wget')]:
    if os.path.exists(wgetPath) == True:
        break
    wgetPath = None

if wgetPath == None:
    print "ERROR: cannot find wget"
    sys.exit()

gunzipPaths = [os.path.join("/","usr","bin","gunzip"),os.path.join("/","bin","gunzip")]
for gunzipPath in gunzipPaths:
    if os.path.exists(gunzipPath) == True:
        break
    gunzipPath = None

if gunzipPath == None:
    print "ERROR: cannot find gunzip exiting..."
    sys.exit()

## functions
def _run_subprocess(cmd):
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

def fetch_file(fetchURL):
    cmd = "%s -N %s"%(wgetPath,fetchURL)
    _run_subprocess(cmd)

def unzip_file(fileName):
    print('unzipping...%s'%fileName)
    cmd = "gunzip -c %s > %s.db"%(fileName,fileName[:-3])
    _run_subprocess(cmd)

def untar_file(fileName):
    print('untarring...%s'%fileName)
    cmd = "tar xzf %s"%fileName
    _run_subprocess(cmd)

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
                "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz",
                "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz",
                uniprotUrl + "idmapping/idmapping.dat.gz",
                uniprotUrl + "idmapping/LICENSE"]

## if blast is enabled then fetch the blast files as well
if CONFIG['blast'] == True:
    filesToFetch.extend([uniprotUrl + "complete/uniprot_sprot.fasta.gz",
                        uniprotUrl + "complete/uniprot_trembl.fasta.gz"])

push_out("fetching blast files = %s"%str(CONFIG['blast']))

for fetchURL in filesToFetch:
    fileName = os.path.split(fetchURL)[-1]
    push_out("fetching %s..."%fileName)
    timeStart = time.time()
    fetch_file(fetchURL)
    fetchTime = time.time() - timeStart
    push_out("...%s"%fetchTime)

    ## unzip the gz files
    if not re.search("\.gz",fileName):
        continue

    if not os.path.exists(fileName[:-3]+".db") or fetchTime > 10:
        if re.search("\.tar\.gz",fileName):
            untar_file(fileName)
        else:
            unzip_file(fileName)

push_out("complete.")
fid.close()

## move back to original directory
os.chdir(cwd)

