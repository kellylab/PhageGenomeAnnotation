#!/bin/bash                                                                                                                                                                
#SBATCH -n 16
#SBATCH -N 1                                                                                                              
#SBATCH -o /home/jbrown/out/171110_all.o      # File to which STDOUT will be written
#SBATCH -e /home/jbrown/out/171110_all.e      # File to which STDERR will be written 
#SBATCH -p sched_mit_chisholm                # name of partition to use  
#SBATCH --mem 60G

module load engaging/anaconda/2.3.0
source activate jb_anaconda
module load engaging/tRNAscan-SE/1.23
module load engaging/prodigal/2.6.3

#python /nobackup1/jbrown/phage_ublast_all.py --outdir /nobackup1/jbrown/newmu/blasts/ --proteindir /nobackup1/jbrown/newmu/proteins/ 1.028.O 1.136.O

python /nobackup1/jbrown/phage_prodigal.py from-genomedir --outdir /nobackup1/jbrown/newmu/ --genomedir /nobackup1/jbrown/newmu/genomes/ 

python /nobackup1/jbrown/phage_ublast_all.py from-proteindir --outdir /nobackup1/jbrown/newmu/blasts/ --proteindir /nobackup1/jbrown/newmu/proteins/ 

python /nobackup1/jbrown/phage_crt.py from-genomedir --outdir /nobackup1/jbrown/newmu/crt/ --genomedir /nobackup1/jbrown/newmu/genomes/ 

python /nobackup1/jbrown/phage_trnascan.py from-genomedir --outdir /nobackup1/jbrown/newmu/trna/ --genomedir /nobackup1/jbrown/newmu/genomes/ 