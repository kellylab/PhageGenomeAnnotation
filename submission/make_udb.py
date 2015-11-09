#!usr/bin/python

import subprocess
import sys

fastafile=sys.argv[1]
databasefile=sys.argv[2]

def run_formatudb(fastafile, databasefile="db.udb", ublast_path="/home/sbiller/usearch7.0.1090_i86linux64"):
    to_run=ublast_path+" -makeudb_ublast "+fastafile+" -output "+databasefile
    subprocess.call(to_run.split(" "))
    
run_formatudb(fastafile, databasefile)