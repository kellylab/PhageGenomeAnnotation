#!/usr/bin/python
import subprocess
import glob
import os
import sys

udb=sys.argv[1]

def run_ublastp(fastafile, out, udb, evalue):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -strand plus -blast6out "+out
    subprocess.call(to_run.split(" "))

phage_prots=glob.glob("/nobackup1/jbrown/annotation/proteins/*.faa")
    
dbname=udb.split("/")[-1].replace(".udb","")

if os.path.exists("/nobackup1/jbrown/annotation/blasts/"+dbname+"/")==False:
    os.mkdir("/nobackup1/jbrown/annotation/blasts/"+dbname+"/")

for p in phage_prots:
    phage=p.split("/")[-1].split("f")[0]
    run_ublastp(fastafile=p, udb=udb, out="/nobackup1/jbrown/annotation/blasts/"+dbname+"/"+phage+"vs."+dbname+".out", evalue="1e-5")
