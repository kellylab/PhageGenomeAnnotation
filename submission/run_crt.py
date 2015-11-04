#!/usr/bin/python

import os
import glob
import subprocess

def run_crt(path_to_crt, input_fasta, output):
    args="java -cp "+path_to_crt+" crt -minNR 2 "+input_fasta+" "+output
    subprocess.call(args.split(" "))

#os.mkdir("./crt")
    
phage_genomes=glob.glob("./genomes/*")

phages=[i.split("/")[-1].replace("final.fasta","") for i in phage_genomes]

path_to_crt="/home/jbrown/programs/CRT1.2-CLI.jar"

for i in range(0, len(phages)):
    input_fasta=phage_genomes[i]
    output="./crt/"+phages[i]+"crt"
    run_crt(path_to_crt, input_fasta, output)
    