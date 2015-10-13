#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -o prod_ublast.out
#SBATCH -e prod_ublast.err
#SBATCH -p sched_mit_hill

module add engaging/python/2.7.9
module add engaging/prodigal/2.6.1

python run_ublasts.py