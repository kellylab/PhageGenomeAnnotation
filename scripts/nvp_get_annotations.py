#!usr/bin/python
import sys
import glob

from nvp_blast_processing import *
from nvp_output_scripts import *

if len(sys.argv)==2:
    phage=sys.argv[1]
else:
    phage="*"
    
phage_genomes=glob.glob("/nobackup1/jbrown/annotation/genomes/"+phage+"final.fasta")

phages=[i.split("/")[-1].replace("final.fasta","") for i in phage_genomes]

big=open("all_annotations.20151109.txt","w")

for p in phages:
    write_gff3_file(p, "/nobackup1/jbrown/annotation/gff3/"+p+"annotations.gff3")
    big.write(cds_blast_annotations_to_table(p))
big.close()
              