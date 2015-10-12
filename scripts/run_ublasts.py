#!/usr/bin/python


import subprocess

def run_ublastp(fastafile, out, udb, evalue):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -blast6out "+out
    subprocess.call(to_run.split(" "))
    
def run_prodigal_phage(inputfasta, out_gene, out_prot):
    to_run="prodigal -i "+inputfasta+" -o "+out_gene+" -a "+out_prot+" -p meta"
    subprocess.call(to_run.split(" "))
    
def run_formatudb(fastafile, databasefile="db.udb"):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -makeudb_ublast "+fastafile+" -output "+databasefile
    subprocess.call(to_run.split(" "))
    
#run_formatudb("./databases/kegg.reduced.fasta","kegg.reduced.udb")

import glob
import os

#udbs=glob.glob("./databases/*.udb")
phage_genomes=glob.glob("./genomes/*")
    
    

#udbs=glob.glob("./databases/*.udb")
#for db in udbs:
 #   dbname=db.split("/")[-1].split(".")[0] 
#     os.mkdir("./blasts/"+dbname)

udbs=["./databases/kegg.reduced.udb"]

for p in phage_genomes:
    phage=p.split("/")[-1].split("f")[0]
    fastafile="./genomes/"+phage+"final.fasta"
    out_gene="./genes/"+phage+"gene"
    out_prot="./proteins/"+phage+"faa"

    #run_prodigal_phage(inputfasta=fastafile, out_gene=out_gene, out_prot=out_prot)

    for db in udbs:
        #dbname=db.split("/")[-1].split(".")[0]              
        #run_ublastp(fastafile=out_prot, udb=db, out="./blasts/"+dbname+"/"+phage+"vs."+dbname+".out", evalue="1e-3")
        dbname="kegg"
        run_ublastp(fastafile=out_prot, udb=db, out="./blasts/"+dbname+"/"+phage+"vs."+dbname+".out", evalue="1e-3")