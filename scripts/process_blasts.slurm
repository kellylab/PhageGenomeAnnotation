#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -o annot_info.out
#SBATCH -e annot_info.err
#SBATCH -p sched_mit_chisholm

module add engaging/python/2.7.9

python get_annotation_info.py