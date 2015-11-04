#!/usr/bin/python


import subprocess
import glob
import os

def run_ublastp(fastafile, out, udb, evalue):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -strand plus -blast6out "+out
    subprocess.call(to_run.split(" "))
    
def run_formatudb(fastafile, databasefile="db.udb"):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -makeudb_ublast "+fastafile+" -output "+databasefile
    subprocess.call(to_run.split(" "))

phage_prots=glob.glob("./proteins/*.faa")
    

os.mkdir("./blasts/tara.translated/")

db="./databases/tara.translated.udb"

for p in phage_prots[0:100]:
    phage=p.split("/")[-1].split("f")[0]
    dbname="tara.translated"
    run_ublastp(fastafile=p, udb=db, out="./blasts/"+dbname+"/"+phage+"vs."+dbname+".out", evalue="1e-5")
