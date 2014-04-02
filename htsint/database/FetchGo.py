#!/usr/bin/env python
"""
fetch gene ontology files
"""
                                                          
import os,sys,subprocess,re,time,csv
from htsint import __basedir__
sys.path.append(__basedir__)

try:
    from configure import CONFIG
except:
    CONFIG = None

if CONFIG == None:
    raise Exception("You must create a configure.py before running FetchGo.py")

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
    print 'unzipping', fileName
    cmd = "gunzip -c %s > %s.db"%(fileName,fileName[:-3])
    _run_subprocess(cmd)

## prepare a log file
fid = open('fetchgo.log','w')
writer = csv.writer(fid)

def push_out(line):
    writer.writerow([line])
    print line

push_out(sys.argv[0])
push_out(time.asctime())
push_out("fetching files...")

## fetch the go term database
filesToFetch = ["go.obo"]
for fileName in filesToFetch:
    urlBase = "http://geneontology.org/ontology/"
    fetchURL = urlBase + fileName
    timeStart = time.time()
    fetch_file(fetchURL)
    fetchTime = time.time() - timeStart

fid.close()

## move back to original directory
os.chdir(cwd)
