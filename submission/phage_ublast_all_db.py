#!/usr/bin/python
import subprocess
import glob
import os
import sys

'''
This script blasts either an individual phage against all annotation databases if included in the script argument, or if no phage is included for the system arguments, the script blasts all phage proteomes against all annotation databases.
'''

if len(sys.argv)==2:
    phage=sys.argv[1]
else:
    phage="*"
    
def run_ublastp(fastafile, out, udb, evalue):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -strand plus -blast6out "+out
    subprocess.call(to_run.split(" "))

phage_prots=glob.glob("/nobackup1/jbrown/annotation/proteins/"+phage+"faa")

udbs=["aclame.udb", "cogs_2003-2014.udb", "CVP.udb", "eggnog4.udb", "Pfam.udb", "pog.udb", "tara.translated.udb","kegg.reduced.fasta"]

dbnames=["aclame","cogs_2003-2004","CVP","eggnog","Pfam","pog","tara.translated","kegg"]

for dbname in dbnames:
    if os.path.exists("/nobackup1/jbrown/annotation/blasts/"+dbname+"/")==False:
        os.mkdir("/nobackup1/jbrown/annotation/blasts/"+dbname+"/")

for p in phage_prots:
    phage=p.split("/")[-1].split("f")[0]
    for i in range(0, len(udbs)):
        udb="/nobackup1/jbrown/annotation/databases/"+udbs[i]
        dbname=dbnames[i]
        run_ublastp(fastafile=p, udb=udb, out="/nobackup1/jbrown/annotation/blasts/"+dbname+"/"+phage+"vs."+dbname+".out", evalue="1e-5")
