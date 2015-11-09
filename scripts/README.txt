## README

#### scripts in this file are for annotation of nahant phage genomes and general functions for manipulation of genomic data.

### Pipeline scripts:

>### nvp_blast_processing.py:

>>set of functions for reading blast files and creating a dictionary of blast annotations from various BLAST comparisons

>### nvp_output_scripts.py:

>>contains functions to write output of BLAST analysis plus additional analyses to a gff3 formatted file, or a general table format.

>### nvp_get_annotations.py:

>>uses functions in the above two files to create annotated phage genome files in .gff3 format as well as ouput a large table of all ORF annotations for the dataset.


>#### run_crt.py:

>>calls and runs crt (found here) which scans genomes for CRISPR-like elements

>#### tRNAscan_nvp_phages.py

>>calls and runs tRNA scan

>#### ublast_phage_proteins.py

>>runs ublast on phage proteins against databases used for annotation

#### Other scripts

>fna_from_prod_and_fasta.py

>>extracts na sequences from a genomic fasta based on standard prodigal output

>split_fasta.py

>>splits fasta file into user-designated number of smaller files