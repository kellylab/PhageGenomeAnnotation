#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -o ublast_egg.out
#SBATCH -e ublast_egg.err
#SBATCH -p sched_mit_hill
#SBATCH --mem 60000

module add engaging/python/2.7.9

python split_fasta.py "./databases/eggnog4.proteins.all.fa" "8"

python eggnog_blasts.py