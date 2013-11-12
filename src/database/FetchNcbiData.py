#!/usr/bin/env python
                                                          
import os,sys,subprocess,re,time,csv

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

def unzip_file(fileName,trim):
    print 'unzipping', fileName
    if trim == None:
        cmd = "tar xzf %s"%(fileName)
    else:
        cmd = "gunzip -c %s > %s.db"%(fileName,fileName[:-trim])
    _run_subprocess(cmd)


## prepare a log file
fid = open('fetchncbi.log','w')
writer = csv.writer(fid)

def push_out(line):
    writer.writerow([line])
    print line

push_out(sys.argv[0])
push_out(time.asctime())
push_out("fetching files...")

## fetch the NCBI databases
filesToFetch = ["taxdump.tar.gz"]

for fileName in filesToFetch:
    urlBase = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/"
    fetchURL = urlBase + fileName
    timeStart = time.time()
    fetch_file(fetchURL)
    fetchTime = time.time() - timeStart

    trim = None
    if fetchTime > 10:
        unzip_file(fileName,trim)

filesToFetch = ["gene2accession.gz",
                "gene_info.gz",
                "gene2go.gz",
                "gene_history.gz"]

for fileName in filesToFetch:
    urlBase = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/"
    fetchURL = urlBase + fileName
    timeStart = time.time()
    fetch_file(fetchURL)
    fetchTime = time.time() - timeStart

    trim = 3
    if  os.path.exists(fileName[:-trim]+".db") == False or fetchTime > 20:
        unzip_file(fileName,trim)

print 'done.'

fid.close()
