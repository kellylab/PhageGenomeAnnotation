#!usr/bin/python
'''
this script runs tRNA scan on all nahant phages, unless one phage is specified.
'''

import subprocess
import glob
import os
import sys


if len(sys.argv)==2:
    phage=sys.argv[1]
else:
    phage="*"
    
if os.path.exists("/nobackup1/jbrown/annotation/trna/")==False:
    os.mkdir("/nobackup1/jbrown/annotation/trna/")

#NOTE: -N flag in tRNA scan outputs codons instead of anti-codons
def run_trna_scan(input_file, output):
    if os.path.exists(output):
        os.remove(output)
    args=["tRNAscan-SE", "-o", output, "-G", "-D", input_file]
    subprocess.call(args)

phage_genomes=glob.glob("/nobackup1/jbrown/annotation/genomes/"+phage+"final.fasta")

phages=[i.split("/")[-1].replace("final.fasta","") for i in phage_genomes]

for p in range(0, len(phage_genomes)):
    input_file=phage_genomes[p]
    phage=phages[p]
    output="/nobackup1/jbrown/annotation/trna/"+phage+"trna"
    run_trna_scan(input_file, output)
    
    