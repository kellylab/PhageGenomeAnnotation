#!usr/bin/python
#The purpose of this script is to grab gene sequences from a genomic contig based on coordinates in prodigal .gff output
#script currently customized for nahant phage genomes... 

from Bio.Seq import reverse_complement
from Bio.Seq import translate
from pyfaidx import Fasta
import re
import sys
'''
prodigal_file=sys.argv[1]
genomic_fasta_file=sys.argv[2]

if len(sys.argv)==4:
    output_file=sys.argv[3]
else:
    output_file=genomic_fasta_file.replace(".fasta",".cds.fna")
'''

class cds(object):
    locus_tag=""
    start=""
    stop=""
    strand=""
    tag=""
    
    def __init__(self, locus_tag, start, stop, strand, tag):
        self.locus_tag=locus_tag
        self.start=start
        self.stop=stop
        self.strand=strand
        self.tag=tag
        
def adjust_plus_strand_position(position):
    position=int(position.replace("<","").replace(">",""))-1
    return position

def adjust_negative_strand_position(position):
    position=int(position.replace("<","").replace(">",""))
    return(position)

def check_partial(start, stop):
    if re.search(r">|<", str(start)+str(stop)):
        tag="partial_cds"
    else:
        tag="cds"
    return tag

def get_digits(prodfile):
    prod=open(prodfile).readlines()
    digits=len(str(len(prod)/2))
    return digits

def set_cds_objects(prodigal_file):      #prod prodigal file

    digits=get_digits(prodigal_file)
    cds_records=[]
    
    prod=open(prodigal_file).readlines()
    Sequence=prod[0].split(";")[2].split("=")[1].replace('"','')
    phage=Sequence.split("_")[1]
    
    for i in range(2,len(prod)-1,2):
        #extract locus tag:
        info=prod[i+1]

        number=info.split(";")[0].split("=")[2].split("_")[1]
        z="0"*(digits-len(number))
        t="NVP"+phage.replace(".","")+"_"+z+number
        
        loc=prod[i]
        if len(loc.split())==2:
            if "complement" in prod[i].split()[1]:
                strand="-"
                start=loc.split("..")[1].replace(")\n","")
                stop=loc.split("(")[1].split("..")[0]
                tag=check_partial(start, stop)
                start=adjust_negative_strand_position(start)
                stop=adjust_negative_strand_position(stop)

            else:
                strand="+"
                start=loc.split()[1].split("..")[0]
                stop=loc.split()[1].split("..")[1].replace("\n","")
                tag=check_partial(start, stop)
                start=adjust_plus_strand_position(start)
                stop=adjust_plus_strand_position(stop)            
        
        obj=cds(t, start, stop, strand, tag)
        cds_records.append(obj)
    return cds_records

def get_na_cds_fasta(genomic_fasta_file, prodigal_file, output_file):
    cds_list=set_cds_objects(prodigal_file)
    gs=Fasta(genomic_fasta_file)
    gseq=str(gs[0])
    
    out=""
    
    for record in cds_list:
        if record.strand=="-":
         
            na_cds_seq=reverse_complement(gseq[record.stop:record.start])
            #print "negative strand translation:" + translate(na_cds_seq)
        else:
            na_cds_seq=gseq[record.start:record.stop]
            #print "positive strand translation:"+translate(na_cds_seq)
            
        out+=">"+record.locus_tag+" "+record.tag+"\n"+na_cds_seq+"\n"
    #return out
    out1=open(output_file, "w")
    out1.write(out)
    out1.close()

