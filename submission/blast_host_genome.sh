#!/bin/bash                                                                                                                                                                
#SBATCH -n 16                                                                                                                                                              
#SBATCH -N 1                                                                                                              
#SBATCH -o 20160316_blast.o      # File to which STDOUT will be written                                                                                   
#SBATCH -e 20160316_blast.e      # File to which STDERR will be written                                                                                   
#SBATCH -p sched_mit_chisholm                # name of partition to use      
#SBATCH --mem 250GB
#SBATCH --exclusive

module load engaging/parallel/20150522

fa=/nobackup1/jbrown/vibrio_genomes/10N.261.54.A9_cds_prod.faa
out=/nobackup1/jbrown/vibrio_genomes/blasts

db1=/nobackup1/jbrown/annotation/databases/Pfam.udb 
db2=/nobackup1/jbrown/annotation/databases/eggnog4.udb  
db3=/nobackup1/jbrown/annotation/databases/kegg_reduced.fasta
 
out1=${out}/10N.261.54.A9_cds_prod_vs_pfam.out
out2=${out}/10N.261.54.A9_cds_prod_vs_eggnog.out
out3=${out}/10N.261.54.A9_cds_prod_vs_kegg.out


/home/sbiller/usearch7.0.1090_i86linux64 -ublast ${fa} -db ${db2} -evalue 1e-5 -accel 0.5 -strand plus -blast6out ${out2}
/home/sbiller/usearch7.0.1090_i86linux64 -ublast ${fa} -db ${db3} -evalue 1e-5 -accel 0.5 -strand plus -blast6out ${out3}
