#!usr/bin/python
import os
import subprocess

def run_prodigal_phage(inputfasta, out_gene):
    to_run="prodigal -i "+inputfasta+" -f gff -o "+out_gene+" -p meta"
    subprocess.call(to_run.split(" "))
    
group15=['2.092.O.', '1.011.O.', '1.141.A.', '1.095.O.', '1.107.B.', '1.069.O.', 
         '1.008.O.', '1.102.O.', '1.107.C.', '1.080.O.', '1.062.O.', '1.043.O.', 
         '1.249.B.', '1.249.A.', '1.020.O.', '1.048.O.', '1.057.O.', '1.125.O.', 
         '1.044.O.', '1.107.A.', '1.040.O.']

#os.mkdir("./10kb_phages/prod_gff/")

for p in group15:
    inputfasta="./10kb_phages/genomes/"+p+"final.fasta"
    out_gene="./10kb_phages/prod_gff/"+p+"prod.gff"
    run_prodigal_phage(inputfasta, out_gene)
    
    