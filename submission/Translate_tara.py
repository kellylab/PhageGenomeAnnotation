#!usr/bin/python

from pyfaidx import Fasta
from Bio.Seq import translate

tara=Fasta("./databases/OM-RGC_seq.release.fna")

tara_aa=open("./databases/OM-RGC_seq.translated.fasta","w")

for s in tara.keys():
    tara_aa.write(">"+s+"\n"+translate(tara[s])+"\n")

tara_aa.close()