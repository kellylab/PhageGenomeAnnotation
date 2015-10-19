#!usr/bin/python

import glob
import subprocess
import os

def run_ublastp(fastafile, out, udb, evalue):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -blast6out "+out
    subprocess.call(to_run.split(" "))

phages=glob.glob("./proteins/*")
eggnog_dbs=glob.glob("./databases/eggnog4.proteins.all.*.fasta")

os.mkdir("./temp/")

for p in phages[1:10]:
    blast_list=["cat"]
    for e in eggnog_dbs:
        fastafile=p
        udb=e
        evalue=1e-5
        out="./temp/"+phage+"vs."+udb.replace("./databases/","").replace(".fasta","")+".out
        run_ublastp(fastafile, out, udb, evalue)
        blast_list.append(out)
    args=blast_list+[">","./blasts/eggnog/"+phage+"vs.eggnog.out"]
    
subprocess.call(["rm","-rf","./temp/"])