#!/usr/bin/python


import subprocess
import glob

def run_ublastp(fastafile, out, udb, evalue):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -strand plus -blast6out "+out
   
    subprocess.call(to_run.split(" "))
    

phage_prots=glob.glob("./proteins/*.faa")
    

db="./databases/tara.translated.udb"

for p in phage_prots[100:200]:
    phage=p.split("/")[-1].split("f")[0]
    dbname="tara.translated"
    run_ublastp(fastafile=p, udb=db, out="./blasts/"+dbname+"/"+phage+"vs."+dbname+".out", evalue="1e-5")
