#!/usr/bin/env python                                                                                                                                                           
import os,sys,subprocess,re,time

baseDir =  os.path.join("..",os.path.realpath(os.path.dirname(__file__)))

## check for the necessary programs
for wgetPath in [os.path.join("/","usr","bin","wget"), os.path.join('/','usr','local','bin','wget')]:
    if os.path.exists(wgetPath) == True:
        break
    wgetPath = None

if wgetPath == None:
    print "ERROR: cannot find wget"
    sys.exit()

for gunzipPath in [os.path.join("/","usr","bin","gunzip"), os.path.join('/','bin','gunzip')]:
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

def unzip_file(fileName,trim):
    filePath = os.path.join(baseDir,fileName)
    print 'unzipping', fileName
    if trim == 3:
        cmd = "%s -c %s > %s.db"%(gunzipPath,filePath,filePath[:-trim])
    else:
        cmd = "%s -c %s > %s.fasta"%(gunzipPath,filePath,filePath[:-trim])
    _run_subprocess(cmd)

## fetch the nr and swissprot databases
filesToFetch = ['nr.gz','swissprot.gz','nt.gz']

for fileName in filesToFetch:
    filePath = os.path.join(baseDir,fileName)
    urlBase = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/"
    fetchURL = urlBase + fileName
    timeStart = time.time()
    fetch_file(fetchURL)
    fetchTime = time.time() - timeStart

    trim = 3
    if  os.path.exists(fileName[:-trim]+".db") == False or fetchTime > 60:
        unzip_file(fileName,trim)

## fetch the uniprot databases
filesToFetch = ["uniprot_trembl.fasta.gz"]

for fileName in filesToFetch:
    filePath = os.path.join(baseDir,fileName)
    urlBase = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
    fetchURL = urlBase + fileName
    timeStart = time.time()
    fetch_file(fetchURL)
    fetchTime = time.time() - timeStart

    trim = 9
    if  os.path.exists(fileName[:-trim]+".fasta") == False or fetchTime > 60:
        unzip_file(fileName,9)

print 'done.'
