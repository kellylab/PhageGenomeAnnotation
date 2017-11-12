#!/bin/bash                                                                                                                                                                
#SBATCH -n 16
#SBATCH -N 1                                                                                                              
#SBATCH -o /home/jbrown/out/171111_all.o      # File to which STDOUT will be written
#SBATCH -e /home/jbrown/out/171111_all.e      # File to which STDERR will be written 
#SBATCH -p sched_mit_chisholm                # name of partition to use  
#SBATCH --mem 60G

module load engaging/anaconda/2.3.0
source activate jb_anaconda

python /nobackup1/jbrown/github/nvp_annotations.py write-gff3s --genome-dir /nobackup1/jbrown/nvp_for_ncbi/genomes/ --outdir /nobackup1/jbrown/nvp_for_ncbi/ \
--ips-dir /nobackup1/jbrown/annotation/protein_ips_calls/ --blast-dir /nobackup1/jbrown/annotation/blasts/ \
--trna-dir /nobackup1/jbrown/annotation/trna/ --crt-dir /nobackup1/jbrown/annotation/crt/ \
--prod-dir /nobackup1/jbrown/annotation/genes/ --prot-dir /nobackup1/jbrown/annotation/proteins/
