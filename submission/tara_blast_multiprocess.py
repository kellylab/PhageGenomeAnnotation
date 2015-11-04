#!usr/bin/python

import glob
import os
from multiprocessing import Pool
import subprocess

def run_formatudb(fastafile, databasefile="db.udb"):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -makeudb_ublast "+fastafile+" -output "+databasefile
    subprocess.call(to_run.split(" "))
    
run_formatudb("./databases/OM-RGC_seq.translated.fasta","./databases/tara.translated.udb")

def run_ublastp(fastafile, out, udb, evalue):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -blast6out "+out
    subprocess.call(to_run.split(" "))

def blast_phage_list(phage):
    phage=p.split("/")[-1].split("f")[0]
    dbname="tara.translated"
    run_ublastp(fastafile=p, udb=db, out="./blasts/"+dbname+"/"+phage+"vs."+dbname+".out", evalue="1e-5")
        

phage_prots=glob.glob("./proteins/*")

os.mkdir("./blasts/tara.translated/")

db="./databases/tara.translated.udb"

#splits=range(0, len(phage_prots), len(phage_prots)/4)
#subset_lists=[]


#for i in range(0, len(splits)-1):
#    subset_list=phage_prots[splits[i]:splits[i+1]]
#    subset_lists.append(subset_list)

pool=Pool(processes=4)
pool.map(blast_phage_list, phage_prots)