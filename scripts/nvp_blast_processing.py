#!usr/bin/python
from __future__ import division
from pyfaidx import Fasta
import re


#paths
prod_path="/nobackup1/jbrown/annotation/genes/"
faa_path="/nobackup1/jbrown/annotation/proteins/"
    
pfam_blast_path="/nobackup1/jbrown/annotation/blasts/Pfam/"
cog_blast_path="/nobackup1/jbrown/annotation/blasts/cogs_2003-2014/"
aclame_blast_path="/nobackup1/jbrown/annotation/blasts/aclame/"
cvp_blast_path="/nobackup1/jbrown/annotation/blasts/CVP/"
kegg_blast_path="/nobackup1/jbrown/annotation/blasts/kegg/"
tara_blast_path="/nobackup1/jbrown/annotation/blasts/tara.translated/"
egg_blast_path="/nobackup1/jbrown/annotation/blasts/eggnog/"
pog_blast_path="/nobackup1/jbrown/annotation/blasts/pog/"


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

#functions for adding annotations/info to BLAST hit based on BLAST database

#functions for adding annotations/info to BLAST hit based on BLAST database
from Bio.KEGG import REST
import sqlite3

###sqlite3 database checking functions:
def query_og1_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT OG from "+tbl+" where ID='"+qid+"'")
    output=c.fetchone()
    result=output[0]
    conn.close()
    return result

def query_og2_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT function from "+tbl+" where OG='"+qid+"'")
    output=c.fetchone()
    result=output[0]
    conn.close()
    return result

def query_func_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT function from "+tbl+" where ID='"+qid+"'")
    output=c.fetchone()
    result=output[0]
    conn.close()
    return result

###find go and functional annotations based on BLAST identifier:
def add_kegg_descript(hit):
    try:
        desc= REST.kegg_find("genes", hit).read()
        K=re.search(r"K[0-9]{5}", desc)
        KEGG=K.group(0)
        a=re.search(r"(?<=K[0-9]{5}).*", desc)
        ann=a.group(0)
        return [KEGG, ann]
    except:
        return ["none", "none"]

def add_cog_descript(hit):
    cog=query_og1_tbl((hit.split("|")[1]),"cog1")
    func=query_og2_tbl(cog,"cog2").replace("\n","")
    return [cog, func]

def add_pfam_descript(hit):
    pfam=query_og1_tbl(hit,"pfam1").split(".")[0]
    function=query_og2_tbl(pfam, "pfam2").replace("\n","")
    return [pfam, function]

def add_aclame_descript(hit):
    annotation=query_func_tbl(hit, "aclame")
    return [hit, annotation]

def add_cvp_descript(hit):
    func=query_func_tbl(hit, "cvp")
    return [hit, func]

##need to add if statement to deal with egg database entries that do not have an OG...
def add_egg_descript(hit):
    try:
        og=query_og1_tbl(hit, "egg1")
        func=query_og2_tbl(og, "egg2")
    except:
        og="none"
        func="none"
    return [og, func]

def add_tara_descript(hit):   #right now just adding the closest hit, TARA sequences come with COG/Pfam info etc 
    return [hit, hit]

db_dict={"kegg":add_kegg_descript, "cog":add_cog_descript, "pfam":add_pfam_descript, "aclame":add_aclame_descript,
        "cvp":add_cvp_descript, "tara":add_tara_descript, "egg":add_egg_descript}

def add_pog_descript(hit):
    og=query_og1_tbl(hit.split("|")[1], "pog")
    function=query_og2_tbl(og, "pog")
    return [og, function]

db_dict={"kegg":add_kegg_descript, "cog":add_cog_descript, "pfam":add_pfam_descript, "aclame":add_aclame_descript,
        "cvp":add_cvp_descript, "tara":add_tara_descript, "egg":add_egg_descript, "pog":add_pog_descript}

def annotated_blast_dict(blast, prod, faa, db, phage):
    blast_dict=set_up_blast_dict(blast, prod, faa, phage)
    blast_db_function=db_dict[db]
    for i in blast_dict.keys():
        hit=blast_dict[i][0]
        info=blast_db_function(hit)
        blast_dict[i]+=info  
    return blast_dict

#load blast files for genome into dict of blast results 
def load_blast_files(phage):
    prod=prod_path+phage+"gene"
    faa=faa_path+phage+"faa"
    pfam_blast=pfam_blast_path+phage+"vs.Pfam.out"
    cog_blast=cog_blast_path+phage+"vs.cogs_2003-2014.out"
    aclame_blast=aclame_blast_path+phage+"vs.aclame.out"
    cvp_blast=cvp_blast_path+phage+"vs.CVP.out"
    kegg_blast=kegg_blast_path+phage+"vs.kegg.out"
    tara_blast=tara_blast_path+phage+"vs.tara.translated.out"
    egg_blast=egg_blast_path+phage+"vs.eggnog.out"
    pog_blast=pog_blast_path+phage+"vs.pog.out"
    
    kegg_blast_dict=annotated_blast_dict(blast=kegg_blast, prod=prod, faa=faa, db="kegg", phage=phage)
    pfam_blast_dict=annotated_blast_dict(blast=pfam_blast, prod=prod, faa=faa, db="pfam", phage=phage)
    cog_blast_dict=annotated_blast_dict(blast=cog_blast, prod=prod, faa=faa, db="cog", phage=phage)
    aclame_blast_dict=annotated_blast_dict(blast=aclame_blast, prod=prod, faa=faa, db="aclame", phage=phage)
    cvp_blast_dict=annotated_blast_dict(blast=cvp_blast, prod=prod, faa=faa, db="cvp", phage=phage)
    tara_blast_dict=annotated_blast_dict(blast=tara_blast, prod=prod, faa=faa, db="tara", phage=phage)
    egg_blast_dict=annotated_blast_dict(blast=egg_blast, prod=prod, faa=faa, db="egg", phage=phage)
    pog_blast_dict=annotated_blast_dict(blast=pog_blast, prod=prod, faa=faa, db="pog", phage=phage)
    
    blasts={"kegg":kegg_blast_dict, "pfam":pfam_blast_dict, "cog":cog_blast_dict, "aclame":aclame_blast_dict,"cvp":cvp_blast_dict, "tara":tara_blast_dict, "egg":egg_blast_dict, "pog":pog_blast_dict}
    return blasts