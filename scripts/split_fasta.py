#!usr/bin/python

import sys
from pyfaidx import Fasta
import numpy as np
import re

#code to split a fasta file into the desired number of fasta files:
fasta_file=sys.argv[1]
number_files=int(sys.argv[2])

def split_fasta(number_files, fasta_file):
    try:
        fasta=Fasta(fasta_file)
    except:
        print "could not open fasta"
        exit()
    number_seqs=len(fasta.keys())
    splits=int(np.ceil(number_seqs/number_files))
    #print(splits)
    ranges=range(0, number_seqs, splits)
    print(ranges)
    ranges[-1]=number_seqs
    print(ranges)

    
    for i in range(0, number_files):
        start=ranges[i]
        stop=ranges[i+1]
        label=re.sub(r"\.fa.*","."+str(i+1)+".fasta", fasta_file)
        out=open(label,"w")
        
        for f in fasta.keys()[start:stop]:
            out.write(">"+f+"\n"+str(fasta[f])+"\n")
        out.close()
        
split_fasta(number_files, fasta_file)