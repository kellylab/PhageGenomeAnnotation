#!usr/bin/python

import os
import shutil

os.mkdir("./10kb_phages/")
os.mkdir("./10kb_phages/genomes/")
os.mkdir("./10kb_phages/prodigal/")
os.mkdir("./10kb_phages/prelim_annotation/")
os.mkdir("./10kb_phages/fna/")
os.mkdir("./10kb_phages/protein/")

group15=['2.092.O.', '1.011.O.', '1.141.A.', '1.095.O.', '1.107.B.', '1.069.O.', 
         '1.008.O.', '1.102.O.', '1.107.C.', '1.080.O.', '1.062.O.', '1.043.O.', 
         '1.249.B.', '1.249.A.', '1.020.O.', '1.048.O.', '1.057.O.', '1.125.O.', 
         '1.044.O.', '1.107.A.', '1.040.O.']

for phage in group15:
    genome="/nobackup1/jbrown/annotation/genomes/"+phage+"final.fasta"
    prod="/nobackup1/jbrown/annotation/genes/"+phage+"gene"
    gff3="/nobackup1/jbrown/annotation/gff3/"+phage+"cds.gff3"
    fna="/nobackup1/jbrown/annotation/fna/"+phage+"fna"
    protein="/nobackup1/jbrown/annotation/proteins/"+phage+"faa"
    
    shutil.copyfile(genome, "./10kb_phages/genomes/"+phage+"final.fasta")
    shutil.copyfile(prod, "./10kb_phages/prodigal/"+phage+"gene")
    shutil.copyfile(fna, "./10kb_phages/fna/"+phage+"fna")
    shutil.copyfile(protein, "./10kb_phages/protein/"+phage+"faa")
    if os.path.exists(gff3):
        shutil.copyfile(gff3, "./10kb_phages/prelim_annotation/"+phage+"cds.gff3")
    
