## README

developed with python 2.7
requires pandas and click

nvp_annotation.py will run entire annotation pipeline

phage_crt.py, phage_ublast.py, phage_prodigal.py and phage_trnascan.py will run specific components

Note:
* Annotation scripts use a sqlite database built on engaging here: /pool001/jbrown/blast_db.sqlite
  * This is currently hard coded into nvp_blast_processing.py
  * Some scripts used to build the sqlite database can be found in sqlite_db_build directory
