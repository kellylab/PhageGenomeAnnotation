#!usr/bin/python

import sqlite3
import os
import glob
from annotation_functions import *

def query_tara_db(tid):
    conn=sqlite3.connect('/pool001/jbrown/tara_db.sqlite')
    c=conn.cursor()

    c.execute("SELECT * from taratbl where ID='"+tid+"'")
    output=c.fetchall()
    gene=output[0][1]
    egg=output[0][2]
    ko=output[0][3]
    kfunc=output[0][4]
    conn.close()
    
    return [gene, egg, ko, kfunc]

phage_genomes=glob.glob("/nobackup1/jbrown/annotation/genomes/*fasta")
phages=[i.split("/")[-1].split("f")[0] for i in phage_genomes]

phage="1.161.O."
 
faa="/nobackup1/jbrown/annotation/proteins/"+phage+"faa"
prod="/nobackup1/jbrown/annotation/genes/"+phage+"gene"

tara="/nobackup1/jbrown/annotation/blasts/tara.translated/"+phage+"vs.tara.translated.out"

out=open("/nobackup1/jbrown/annotation/tara_annotations/"+phage+"tara.txt","w")

tarablast=set_up_blast_dict(tara, prod, faa, phage)

for k in tarablast.keys():
    tid=tarablast[k][0]
    tann=query_tara_db(tid)
    together=tarablast[k]+tann
    together=map(str, together)
    to_write="\t".join(together)
    out.write(k+"\t"+to_write+"\n")
out.close()
              
    

    
