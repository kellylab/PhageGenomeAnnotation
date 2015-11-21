#!usr/bin/python
import subprocess


def run_ublastp(fastafile, out_file, udb, evalue, ublast_path="/home/sbiller/usearch7.0.1090_i86linux64"):
    to_run=ublast_path+" -ublast "+fastafile+" -db "+udb+" -evalue "+str(evalue)+" -accel 0.5 -blast6out "+out_file
    subprocess.call(to_run.split(" "))
    
run_ublastp("allphage_goodProteins.fasta", "phage_vs_host.out","allhost_goodProteins.fasta", 1e-3)