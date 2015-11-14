#!usr/bin/python
import sys
import glob

from nvp_blast_processing import *
from nvp_output_scripts import *

output=sys.argv[1]

if len(sys.argv)==3:
    phage=sys.argv[2]
else:
    phage="*"
    
phage_genomes=glob.glob("/nobackup1/jbrown/annotation/genomes/"+phage+"final.fasta")

phages=[i.split("/")[-1].replace("final.fasta","") for i in phage_genomes]

out=open(output, "w")

for p in phages:
    out.write(kegg_egg_pog_tbl(p))
    
out.close()