#!usr/bin/python

import glob
import os
from fna_from_prod_and_fasta import *

#os.mkdir("./fna")

#phage_genomes=glob.glob("./genomes/*")

#phages=[i.split("/")[2].split("f")[0] for i in phage_genomes]

phages="1.141.A."

for phage in phages:
    prod="./genes/"+phage+"gene"
    genome="./genomes/"+phage+"final.fasta"
    outfile="./fna/"+phage+"fna"
    
    get_na_cds_fasta(genome, prod, outfile)
    