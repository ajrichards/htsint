#!/bin/bash
#$ -S /bin/bash
#$ -j yes
#$ -o /usr/share/blast/fetch.log

# My command lines I want to run on the cluster
/usr/bin/python /usr/share/blast/fetchBlastDBs.py