#Functions to run annotation-associated scripts:

import subprocess
import os
import cPickle as pickle
from __future__ import division
from Bio.KEGG import REST
from pyfaidx import Fasta
import re
import sqlite3

#paths
prod_path="/nobackup1/jbrown/annotation/gene/"
faa_path="/nobackup1/jbrown/annotation/proteins/"
    
pfam_blast_path="/nobackup1/jbrown/annotation/blasts/Pfam/"
cog_blast_path="/nobackup1/jbrown/annotation/blasts/cogs_2003-2014/"
aclame_blast_path="/nobackup1/jbrown/annotation/blasts/aclame/"
cvp_blast_path="/nobackup1/jbrown/annotation/blasts/CVP/"
kegg_blast_path="/nobackup1/jbrown/annotation/blasts/kegg/"
tara_blast_path="/nobackup1/jbrown/annotation/blasts/tara.translated/"


def run_prodigal_phage(inputfasta, out_gene, out_prot):
    to_run="prodigal -i "+inputfasta+" -o "+out_gene+" -a "+out_prot+" -p meta"
    subprocess.call(to_run.split(" "))
    
def run_formatudb(fastafile, databasefile="db.udb", ublast_path="/home/sbiller/usearch7.0.1090_i86linux64"):
    to_run=ublast_path+" -makeudb_ublast "+fastafile+" -output "+databasefile
    subprocess.call(to_run.split(" "))
    
def run_ublastp(fastafile, out_file, udb, evalue, ublast_path="/home/sbiller/usearch7.0.1090_i86linux64"):
    to_run=ublast_path+" -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -blast6out "+out_file+" -top_hit_only"
    subprocess.call(to_run.split(" "))
    
def run_trna_scan(input_file, output):
    if os.path.exists(output):
        os.remove(output)
    args=["tRNAscan-SE", "-o", output, "-G", "-D","-N", input_file]
    subprocess.call(args)
    print("tRNA scan of "+input_file+" is done!")

def run_crt(path_to_crt, input_fasta, output):
    args="java -cp "+path_to_crt+" crt -minNR 2 "+input_fasta+" "+output
    subprocess.call(args.split(" "))
    
#general info

#get number of genes in genome to know how many digits to use in locus tag


def get_digits(faa):
    faa=open(faa).read()
    digits=len(str(faa.count(">")))
    return digits

#create locus tag from protein sequence name in BLAST output file:
def get_locus_tag(line, digits, phage):
    query=line.split("\t")[0].split(" ")[0]
    number=query.split("_")[-1]
    z="0"*(digits-len(number))
    return "NVP"+phage.replace(".","")+"_"+z+number
        

from pyfaidx import Fasta

def get_prot_lens(faa_file, phage):
    len_dict={}
    digits=get_digits(faa_file)
    #def make_seq_len_dict(faa):
    f=Fasta(faa_file)
    for i in f.keys():
        name=get_locus_tag(i, digits=digits, phage=phage)
        length=len(str(f[i]))
        len_dict[name]=length
    return len_dict

#set up dict of general info from BLAST:
def set_up_blast_dict(blast, prod, faa, phage):
    digits=get_digits(faa)
    len_dict=get_prot_lens(faa, phage)
    records=[]
    blast_dict={}
    
    blast=open(blast).readlines()
    for line in blast:
        name=line.split(" ")[0]
        hit=line.split("\t")[1]
        lt=get_locus_tag(name, digits=digits, phage=phage)
        prot_len=len_dict[lt]
        aln_len=int(line.split("\t")[3])
        pct_id=float(line.split("\t")[2])
        ev=line.split("\t")[-2]
        pct_cov=(aln_len/prot_len)*100

        if pct_id>35 and pct_cov>75 and lt not in records:
            records.append(lt)
            blast_dict[lt]=[hit, pct_cov, pct_id, ev]
    return blast_dict

#load blast database dictionaries:
if os.path.exists("/nobackup1/jbrown/annotation/databases/pickled_dicts/"):
    #path to dictionaries stored on server:
    aclame_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/aclame_dict.p","rb"))
    cog_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/cog_dict.p","rb"))
    cog_defs=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/cog_def.p","rb"))
    pfam_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/pfam_dict.p","rb"))
    pfam_defs=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/pfam_def.p","rb"))
    cvp_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/cvp_dict.p","rb"))
    dict_check=True
elif os.path.exists("./databases/pickled_dicts/"):
    #path to dictionaries stored on jmb@alarmism.einstein.yu.edu
    aclame_dict=pickle.load(open("./databases/pickled_dicts/aclame_dict.p","rb"))
    cog_dict=pickle.load(open("./databases/pickled_dicts/cog_dict.p","rb"))
    cog_defs=pickle.load(open("./databases/pickled_dicts/cog_def.p","rb"))
    pfam_dict=pickle.load(open("./databases/pickled_dicts/pfam_dict.p","rb"))
    pfam_defs=pickle.load(open("./databases/pickled_dicts/pfam_def.p","rb"))
    cvp_dict=pickle.load(open("./databases/pickled_dicts/cvp_dict.p","rb"))
    dict_check=True
else:
    print("blast database dictionaries not found")
    dict_check=False

#functions for adding annotations/info to BLAST hit based on BLAST database

def add_kegg_descript(hit):
    
    desc= REST.kegg_find("genes", hit).read()
    K=re.search(r"K[0-9]{5}", desc)
    KEGG=K.group(0)
    a=re.search(r"(?<=K[0-9]{5}).*", desc)
    ann=a.group(0)
    return [KEGG, ann]

def add_cog_descript(hit):
    if dict_check:
        cog=cog_dict[(hit.split("|")[1])]
        func=cog_defs[cog].replace("\n","")
        return [cog, func]
    else:
        print "COG database is not loaded"
        return ""

def add_pfam_descript(hit):
    if dict_check:
        pfam=pfam_dict[hit].split(".")[0]
        function=pfam_defs[pfam].replace("\n","")
        return [pfam, function]
    else:
        print "pfam database is not loaded"
        return ""

def add_aclame_descript(hit):
    if dict_check:
        annotation=aclame_dict[hit]
        return [hit, annotation]
    else:
        print "aclame database not loaded"
        return ""

def add_cvp_descript(hit):
    func=cvp_dict[hit]
    return [hit, func]

def add_tara_descript(hit):   #right now just adding the closest hit, TARA sequences come with COG/Pfam info etc 
    return [hit, hit]

db_dict={"kegg":add_kegg_descript, "cog":add_cog_descript, "pfam":add_pfam_descript, "aclame":add_aclame_descript,
        "cvp":add_cvp_descript, "tara":add_tara_descript}

def annotated_blast_dict(blast, prod, faa, db, phage):
    blast_dict=set_up_blast_dict(blast, prod, faa, phage)
    blast_db_function=db_dict[db]
    for i in blast_dict.keys():
        hit=blast_dict[i][0]
        info=blast_db_function(hit)
        blast_dict[i]+=info
    
    return blast_dict

#i=prodigal line that begins with a location identifier..
#function meant to iterate over the length of a prodigal file every two lines starting at line 3 as such: 
'''
for i in range (2, len(open(prod_file).readlines())-1,2):
    get_prod_cds_info(i,...)
'''

def get_prod_cds_info(i, prod, digits, phage):  
    loc=prod[i]
    if len(loc.split())==2:
        if "complement" in prod[i].split()[1]:
            strand="-"
            start=loc.split("..")[1].replace(")\n","")
            stop=loc.split("(")[1].split("..")[0]
        else:
            strand="+"
            start=loc.split()[1].split("..")[0]
            stop=loc.split()[1].split("..")[1].replace("\n","")
        start=start.replace(">","").replace("<","")
        stop=stop.replace(">","").replace("<","")
    info=prod[i+1]
    number=info.split(";")[0].split("=")[2].split("_")[1]
    z="0"*(digits-len(number))
    t="NVP"+phage.replace(".","")+"_"+z+number
    return [t, start, stop, strand]

#below: considers hits to more informative databases before less informative databases
#dict_list* are lists of blast_dicts and dl*_names are the names of the dicts in the same order

def find_best_hit2(gene_id, dict_list1, dl1_names, dict_list2=[], dl2_names=[]):
    evals=1
    annotation=""
    best_hit=""
    hits=[]
    es=[]
    names=[]
    for i in range(0, len(dict_list1)):
        if gene_id in dict_list1[i].keys():
            hit=dict_list1[i][gene_id]
            hits.append(hit[-1])
            es.append(float(hit[3]))
            names.append(dl1_names[i])
    if len(hits)>0:
        best_annotation=[hits[es.index(min(es))],names[es.index(min(es))]]
    else:
        for i in range(0, len(dict_list2)):
            if gene_id in dict_list2[i].keys():
                hit=dict_list2[i][gene_id]
                hits.append(hit[-1])
                es.append(float(hit[3]))
                names.append(dl2_names[i])
        if len(hits)>0:
            best_annotation=[hits[es.index(min(es))],names[es.index(min(es))]]
        else:
            best_annotation=["",""]

    #print("best annotation for"+gene_id+" is from "+best_hit+" with e-value "+str(evals)+" and annotation of "+annotation)
    return best_annotation 

def gff3_header(prod):
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    return Sequence+"\n"

#merge BLAST results into one gff3
    
def cds_blast_annotations_to_gff3(phage):
    #prodigal and fasta files:
    prod=prod_path+phage+"gene"
    faa=faa_path+phage+"faa"
    
    #blast files:
    pfam_blast=pfam_blast_path+phage+"vs.Pfam.out"
    cog_blast=cog_blast_path+phage+"vs.cogs_2003-2014.out"
    aclame_blast=aclame_blast_path+phage+"vs.aclame.out"
    cvp_blast=cvp_blast_path+phage+"vs.CVP.out"
    kegg_blast=kegg_blast_path+phage+"vs.kegg.out"
    tara_blast=tara_blast_path+phage+"vs.tara.translated.out"

    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    
    #set up dicts from all BLASTs
    kegg_blast_dict=annotated_blast_dict(blast=kegg_blast, prod=prod, faa=faa, db="kegg", phage=phage)
    pfam_blast_dict=annotated_blast_dict(blast=pfam_blast, prod=prod, faa=faa, db="pfam", phage=phage)
    cog_blast_dict=annotated_blast_dict(blast=cog_blast, prod=prod, faa=faa, db="cog", phage=phage)
    aclame_blast_dict=annotated_blast_dict(blast=aclame_blast, prod=prod, faa=faa, db="aclame", phage=phage)
    cvp_blast_dict=annotated_blast_dict(blast=cvp_blast, prod=prod, faa=faa, db="cvp", phage=phage)
    tara_blast_dict=annotated_blast_dict(blast=tara_blast, prod=prod, faa=faa, db="tara", phage=phage)

    #prioritize and name dicts:
    #preferred blast dbs to annotate from if there's a match:
    first_looks=[kegg_blast_dict, pfam_blast_dict, cog_blast_dict, aclame_blast_dict]
    flnames=["kegg","pfam","cog","aclame"]
    #secondary database(s) to annotate from:
    second_look=[cvp_blast_dict]
    slnames=["CVP"]
    #databases with orthologous groups to include in annotation
    OGs=[kegg_blast_dict, pfam_blast_dict, cog_blast_dict]
    OG_names=["KEGG","PFam","COG"]
    #databases where the closest hit will be referenced, but no other info will be provided:
    annotes=[aclame_blast_dict, cvp_blast_dict, tara_blast_dict]
    annotes_names=["ACLAME","CAMERA_viral_proteins","TARA_Oceans_Dataset"]
    
    out=""  #set up string to write to
    
    #run through annotations of each prodigal-identified CDS:
    prod=open(prod).readlines()
    digits=get_digits(faa)
    
    #write gff3 lines from prodigal files and blast dicts:
    for i in range(2,len(prod)-1,2):

        coords=get_prod_cds_info(i, prod, digits, phage)
        locus_tag=coords[0]
        start=coords[1]
        stop=coords[2]
        strand=coords[3]
        
        #set up col9
        col9="ID="+locus_tag
        
        #ID best hit:
        best_hits=find_best_hit2(locus_tag, dict_list1=first_looks, dl1_names=flnames, dict_list2=second_look, dl2_names=slnames)

        #establish name:
        Name=best_hits[0]
        if len(Name)==0:
            col9+='; Name=hypothetical protein'
        else:
            col9+="; Name="+Name.replace('"','')

        #Add OG annotations:
        for d in range(0, len(OGs)):
            og_dict=OGs[d]
            if locus_tag in og_dict.keys():
                col9+='; Ontology_term="'+OG_names[d]+":"+og_dict[locus_tag][-2]+'"'

        #Add db closest hits to notes
        for d in range(0, len(annotes)):
            annote_dict=annotes[d]
            if locus_tag in annote_dict.keys():
                col9+='; note="'+annotes_names[d]+"_best_match:"+annote_dict[locus_tag][-2]+'"'
        out+=Sequence+"\t"+"prod"+"\t"+"CDS"+"\t"+start+"\t"+stop+"\t"+"."+"\t"+strand+"\t"+"0"+"\t"+col9+"\n"

    return out

def CRISPR_gff3(input_fasta, crt_output):
    crtout=open(crt_output).readlines()

    sequence=open(input_fasta).readlines()
    name= [i.replace(">","") for i in sequence if i.startswith(">")][0]
    out=""
    for line in crtout:
        if line.startswith("CRISPR"):
            vec=line.split()
            number=vec[1]
            start=vec[3]
            stop=vec[5]
            ID="NVP"+name.split("_")[1].replace(".","")+"_CRISPR-like_"+number
            out+=name.replace("\n","")+"\t"+"crt"+"\tputative CRISPR feature\t%s\t%s\t.\t.\t.\tID=%s" % (start, stop, ID)
            out+=", note=CRISPR region\n"
    return out

def tRNA_scan_to_gff3(tRNAScanSE_file):
    if os.path.getsize(tRNAScanSE_file)>0:
        t=open(tRNAScanSE_file).readlines()
        tanns=""
        for line in t[3:]:
            l=line.split("\t")

            locus_tag="NVP"+l[0].split("_")[1].replace(".","")+"_tRNA_"+l[1]
            start=l[2]
            stop=l[3]
            if start<stop:
                strand="+"
            else:
                strand="-"
            aa=l[4]
            codon=l[5]
            SeqID=l[0]
            col9="ID="+locus_tag+", aa="+aa+", codon="+codon
            out=SeqID+"\t"+"tRNAScanSE"+"\t"+"tRNA"+"\t"+start+"\t"+stop+"\t"+l[-1].replace("\n","")+"\t"+strand+"\t"+"0"+"\t"+col9+"\n"
            tanns+=out
        return tanns
    else:
        print "no tRNAs found in genome"
        return ""
    
#put them all together:
def write_gff3_file(phage):
    prod=prod_path+phage+".gene"
    faa=faa_path+phage+".faa"
    genomic_fasta="./genomes/%sfinal.fasta" % phage
    
    out=open(phage+"test.gff3","w")
    #out.write(gff3_header(phage+"gene"))
    out.write(cds_blast_annotations_to_gff3(phage))
    
    trna="../tRNA_info/data/nahant_tRNA_count/%strnas.txt" % phage 
    if os.path.getsize(trna)>0:
        out.write(tRNA_scan_to_gff3(trna))
    
    crt_output=phage+"crt"
    out.write(CRISPR_gff3(genomic_fasta, crt_output))
    
    out.close()

def query_tara_db(tid):
    conn=sqlite3.connect('/pool001/jbrown/tara_db.sqlite')
    c=conn.cursor()

    c.execute("SELECT * from taratbl where ID='"+tid+"'")
    output=c.fetchall()
    gene=output[0][1]
    egg=output[0][2]
    ko=output[0][3]
    kfunc=output[0][4]
    conn.close()
    
    return [gene, egg, ko, kfunc]

def cds_blast_annotations_to_table(phage):
    #set up reference files:
    prod=prod_path+phage+"gene"
    faa=faa_path+phage+"faa"
    
    #blast files:
    pfam_blast=pfam_blast_path+phage+"vs.Pfam.out"
    cog_blast=cog_blast_path+phage+"vs.cogs_2003-2014.out"
    aclame_blast=aclame_blast_path+phage+"vs.aclame.out"
    cvp_blast=cvp_blast_path+phage+"vs.CVP.out"
    kegg_blast=kegg_blast_path+phage+"vs.kegg.out"
    tara_blast=tara_blast_path+phage+"vs.tara.translated.out"

    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    
    #set up dicts from all BLASTs
    kegg_blast_dict=annotated_blast_dict(blast=kegg_blast, prod=prod, faa=faa, db="kegg", phage=phage)
    pfam_blast_dict=annotated_blast_dict(blast=pfam_blast, prod=prod, faa=faa, db="pfam", phage=phage)
    cog_blast_dict=annotated_blast_dict(blast=cog_blast, prod=prod, faa=faa, db="cog", phage=phage)
    aclame_blast_dict=annotated_blast_dict(blast=aclame_blast, prod=prod, faa=faa, db="aclame", phage=phage)
    cvp_blast_dict=annotated_blast_dict(blast=cvp_blast, prod=prod, faa=faa, db="cvp", phage=phage)
    tara_blast_dict=annotated_blast_dict(blast=tara_blast, prod=prod, faa=faa, db="tara", phage=phage)

    #prioritize and name dicts:
    #preferred blast dbs to annotate from if there's a match:
    first_looks=[kegg_blast_dict, pfam_blast_dict, cog_blast_dict, aclame_blast_dict]
    flnames=["kegg","pfam","cog","aclame"]
    #secondary database(s) to annotate from:
    second_look=[cvp_blast_dict]
    slnames=["CVP"]
    #databases with orthologous groups to include in annotation
    OGs=[kegg_blast_dict, pfam_blast_dict, cog_blast_dict]
    OG_names=["KEGG","PFam","COG"]
    #databases where the closest hit will be referenced, but no other info will be provided:
    annotes=[aclame_blast_dict, cvp_blast_dict, tara_blast_dict]
    annotes_names=["ACLAME","CAMERA_viral_proteins","TARA_Oceans_Dataset"]
    
    out=open(phage+"annotations.txt","w")
    out.write("Genome\tlocus_tag\ttype\tstart\tstop\tstrand\tbest_hit_annotation\t")
    out.write("KEGG\tPFam\tCOG\tACLAME\tCVP\tTARA\n")
    
    #run through annotations of each prodigal-identified CDS:
    prod=open(prod).readlines()
    digits=get_digits(faa)
    
    #write gff3 lines from prodigal files and blast dicts:
    for i in range(2,len(prod)-1,2):
        out.write(Sequence+"\t")
        
        
        coords=get_prod_cds_info(i, prod, digits, phage)
        locus_tag=coords[0]
        start=coords[1]
        stop=coords[2]
        strand=coords[3]
        
        out.write(locus_tag+"\tcds\t"+start+"\t"+stop+"\t"+strand+"\t")
        
        #ID best hit:
        best_hits=find_best_hit2(locus_tag, dict_list1=first_looks, dl1_names=flnames, dict_list2=second_look, dl2_names=slnames)

        #establish name:
        Annotation=best_hits[0]
        if len(Annotation)==0:
            Annotation='hypothetical protein'
        else:
            Annotation=Annotation.replace('"','')
        out.write(Annotation+"\t")
        
        #Add OG annotations:
        for d in range(0, len(OGs)):
            og_dict=OGs[d]
            if locus_tag in og_dict.keys():
                out.write(og_dict[locus_tag][-2]+"\t")
            else:
                out.write(" \t")

        #Add db closest hits to notes
        for d in range(0, len(annotes)):
            annote_dict=annotes[d]
            if locus_tag in annote_dict.keys():
                out.write(annote_dict[locus_tag][-2]+'\t')
            else:
                out.write(" \t")
        out.write("\n")

    out.close()