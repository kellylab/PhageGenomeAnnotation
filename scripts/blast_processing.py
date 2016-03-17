#!usr/bin/env python
'''
ultimately creates table of annotations 
'''
from __future__ import division
from __future__ import print_function
from pyfaidx import Fasta
import re
from Bio.KEGG import REST
import sqlite3
import itertools
import os.path as op
import subprocess
#import click

pfam_udb = "/nobackup1/jbrown/annotation/databases/Pfam.udb"
eggnog_udb = "/nobackup1/jbrown/annotation/databases/eggnog4.udb"
kegg_udb = "/nobackup1/jbrown/annotation/databases/kegg_reduced.fasta"

def run_ublastp(fastafile, out, udb, evalue=1e-5):
    to_run="/home/sbiller/usearch7.0.1090_i86linux64 -ublast %s -db %s -evalue %s -accel 0.5 -strand plus -blast6out %s" % (fastafile, udb, evalue, out)
    subprocess.call(to_run.split(" "))

udbs = [kegg_udb, eggnog_udb, pfam_udb]

def ublast_udbs(faa, outdir, udbs=[kegg_udb, eggnog_udb, pfam_udb]):
    blast_files = []
    for d in udbs:
        out = op.join(outdir, "%s_vs_%s.out" % (op.basename(faa).split(".")[0], op.basename(d).split(".")[0]))
        if op.exists(out):
            blast_files.append(out)
        else:
        	run_ublastp(faa, out, d)
        	blast_files.append(out)
    return blast_files

def read_fasta(file_handle):
    '''fasta generator'''
    for header, group in itertools.groupby(file_handle, lambda line: line[0] == '>'):
        if header:
            line = group.next()
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq
            
def get_prot_lens(faa_file, phage):
    len_dict = {}
    with open(faa_file) as infile:
        for name, seq in read_fasta(infile):
            len_dict[name] = len(seq)
    return len_dict
    
#set up dict of general info from BLAST:
def set_up_blast_dict(blast, faa, cov_thresh = 75):
    '''set up dict of general info from blast
    Args:
        blast (path): path to blast file
        cov_thresh (numeric): coverage threshold for which to keep a hit.  Default is 75% id
    Returns:
        blast dictionary with sequence name as key, value of [hit, pct_coverage, pct_id, evalue]
    '''

    len_dict = get_prot_lens(faa, phage)
    records = []
    blast_dict = {}
    try:
        with open(blast) as infile:
			for line in infile:
				name = line.split(" ")[0]
				hit = line.split("\t")[1]	
				lt = name    
				prot_len = len_dict[lt]
				aln_len = int(line.split("\t")[3])
				pct_id = float(line.split("\t")[2])
				ev = line.split("\t")[-2]
				pct_cov = (aln_len/prot_len)*100
				if pct_id > 35 and pct_cov > cov_thresh and lt not in records:
					records.append(lt)
					blast_dict[lt] = [hit, pct_cov, pct_id, ev]
    except:
        print "could not open %s" %blast
    return blast_dict


#query sqlite3 table queries for various outputs:

def query_og1_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    '''
    input query id, get orthologous group id
    '''
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT OG from "+tbl+" where ID='"+qid+"'")
    output=c.fetchone()
    result=output[0]
    conn.close()
    return result

def query_og2_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    '''input orthologous group, return function'''
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT function from "+tbl+" where OG='"+qid+"'")
    output=c.fetchone()
    result=output[0]
    conn.close()
    return result

def query_og3_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    '''input id, return category'''
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT category from "+tbl+" where ID='"+qid+"'")
    output=c.fetchone()
    result=output[0]
    conn.close()
    return result

def query_func_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    '''input id, return function'''
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT function from "+tbl+" where ID='"+qid+"'")
    output=c.fetchone()
    result=output[0]
    conn.close()
    return result

def query_phy_tbl(qid, tbl, db_location='/pool001/jbrown/blast_db.sqlite'):
    '''input id, return phylogeny list'''
    conn=sqlite3.connect(db_location)
    c=conn.cursor()
    c.execute("SELECT phy1, phy2, phy3 from "+tbl+" where ID='"+qid+"'")
    output=c.fetchone()
    result=output
    conn.close()
    return result


def reduce_func_len(function):
    '''shorten excessively long protein descriptions'''
    if len(function)>50:
        if function.split(",")>1:
            function=function.split(",")[0]
        elif function.split(";")>1:
            function=function.split(";")[0]
        elif function.split(".")>1:
            function=function.split(".")[0]
        else:
            function=function[0:50]
    else:
        function=function
    return function

###find go and functional annotations based on BLAST identifier:
def add_kegg_descript(hit):
    try:
        desc= REST.kegg_find("genes", hit).read()
        try:
            K=re.search(r"K[0-9]{5}", desc)
            KEGG=K.group(0)
        except:
            KEGG="none"
        try:
            a=re.search(r"(?<=K[0-9]{5}).*", desc).replace("\n","")
            ann=a.group(0)
        except:
            try:
                ann=desc.split("\t")[1].split(";")[0].replace("\n","")
            except:
                ann="none"
    except:
        ann="none"
        KEGG="none"
    ann=reduce_func_len(ann)
    return strip_lines_list([KEGG, ann])

# this is definitely the slowest step of the pipeline    
def add_kegg_descript2(hit):
    '''use KEGG REST api to extract KEGG definitions'''
    try:
        desc= REST.kegg_find("genes", hit).read()
        try:
            K=re.search(r"K[0-9]{5}", desc)
            KEGG=K.group(0)
        except:
            KEGG="none"
        try:
            a=re.search(r"(?<=K[0-9]{5}).*", desc).replace("\n","")
            ann=a.group(0)
        except:
            try:
                ann=desc.split("\t")[1].split(";")[0].replace("\n","")
            except:
                ann="none"
        try:
            mod=REST.kegg_link('module', hit).read()
            module=mod.split(":")[2].split("_")[-1].replace("\n","")
        except:
            module="none" 
    except:
        module="none"
        KEGG="none"
        ann="none"
    ann=reduce_func_len(ann)
    return strip_lines_list([module, KEGG, ann])

def add_cog_descript(hit):
    cog=query_og1_tbl((hit.split("|")[1]),"cog1")
    func=query_og2_tbl(cog,"cog2").replace("\n","")
    func=reduce_func_len(func)
    return strip_lines_list([cog, func])

def add_pfam_descript(hit):
    pfam=query_og1_tbl(hit,"pfam1").split(".")[0]
    function=query_og2_tbl(pfam, "pfam2").replace("\n","")
    function=reduce_func_len(function)
    return strip_lines_list([pfam, function])

def add_aclame_descript(hit):
    annotation=query_func_tbl(hit, "aclame")
    annotation=reduce_func_len(annotation)
    return strip_lines_list([hit, annotation])

def add_cvp_descript(hit):
    func=query_func_tbl(hit, "cvp")
    func=reduce_func_len(func)
    return strip_lines_list([hit, func])

def add_egg_descript(hit):
    try:
        og=query_og1_tbl(hit, "egg1")
        func=query_og2_tbl(og, "egg2")
    except:
        og="none"
        func="none"
    func=reduce_func_len(func)
    return strip_lines_list([og, func])

def add_egg_descript2(hit):
    try:
        og=query_og1_tbl(hit, "egg1")
        func=query_og2_tbl(og, "egg2")
        cat=query_og3_tbl(hit, "egg1")
    except:
        og="none"
        func="none"
        cat="none"
    func=reduce_func_len(func)
    return strip_lines_list([cat, og, func])

def add_tara_descript(hit):   #right now just adding the closest hit, TARA sequences come with COG/Pfam info etc 
    return strip_lines_list([hit, hit])

def add_pog_descript(hit):
    og=query_og1_tbl(hit.split("|")[1], "pog")
    function=query_og2_tbl(og, "pog")
    function=reduce_func_len(function)
    return strip_lines_list([og, function])

def add_pog_descript2(hit):
    og=query_og1_tbl(hit.split("|")[1], "pog")
    function=query_og2_tbl(og, "pog")
    function=reduce_func_len(function)
    phylog=query_phy_tbl(hit.split("|")[1],"pog")
    phylist=[i.split("]")[-1] for i in phylog]
    return strip_lines_list([og, function]+phylist)
    
def annotated_blast_dict(blast, prod, faa, db, phage, cov_thresh=75):
    db_dict={"kegg":add_kegg_descript, 
         "cog":add_cog_descript, 
         "pfam":add_pfam_descript, 
         "aclame":add_aclame_descript,
         "cvp":add_cvp_descript, 
         "tara":add_tara_descript, 
         "egg":add_egg_descript, 
         "pog":add_pog_descript}

    blast_dict=set_up_blast_dict(blast, prod, faa, phage, cov_thresh=cov_thresh)
    blast_db_function=db_dict[db]
    for i in blast_dict.keys():
        hit=blast_dict[i][0]
        info=blast_db_function(hit)
        blast_dict[i]+=info  
    return blast_dict

def enhanced_blast_dict(blast, faa, db, cov_thresh=75):
    db_dict={"kegg":add_kegg_descript2, 
         "egg":add_egg_descript2, 
         "pfam":add_pfam_descript}

    blast_dict=set_up_blast_dict(blast, faa, cov_thresh=cov_thresh)
    blast_db_function=db_dict[db]
    for i in blast_dict.keys():
        hit=blast_dict[i][0]
        info=blast_db_function(hit)
        blast_dict[i]+=info  
    return blast_dict

    
def load_kegg_egg_pfam_blast(blast_list, faa, cov_thresh=75):
    kegg_blast = blast_list[0]
    egg_blast = blast_list[1]
    pfam_blast = blast_list[2]
    kegg_blast_dict = enhanced_blast_dict(blast=kegg_blast, faa=faa, db="kegg", cov_thresh=cov_thresh)
    egg_blast_dict = enhanced_blast_dict(blast=egg_blast, faa=faa, db="egg", cov_thresh=cov_thresh)
    pog_blast_dict = enhanced_blast_dict(blast=pfam_blast, faa=faa, db="pfam", cov_thresh=cov_thresh)
    blasts={"kegg":kegg_blast_dict, 
            "egg":egg_blast_dict, 
            "pfam":pfam_blast_dict}
    return blasts
    

def kegg_egg_pfam_tbl(prefix, outfile, cov_thresh=75):
    #prodigal and fasta files:
    prod=prod_path+phage+"gene"
    faa=faa_path+phage+"faa"
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    
    ogs=["kegg","egg","pfam"]
    og_lens={"kegg":7,"egg":7,"pog":6}
    blast_dict=load_kegg_egg_pfam_blast(phage, cov_thresh=cov_thresh)
    
    with open(outfile, "w") as oh:
        for i in blast_dict.keys():
    #write lines from blast dicts:
    for i in range(2,len(prod)-1,2):
        #col1
        out+=Sequence+"\t"
        
        #Add OG annotations:
        for d in range(0, len(ogs)):
            og_dict=blast_dict[ogs[d]]
            if locus_tag in og_dict.keys():
                for entry in og_dict[locus_tag]:
                    out+=str(entry)+"\t"
            else:
                out+="NA\t"*og_lens[ogs[d]]
        out+="\n"

    return out.replace("\t\n","\n")