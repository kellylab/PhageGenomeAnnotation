#!/usr/bin/python
#first go at script to write a gff3 file output from BLAST results and prodigal calls... better version in development.

import subprocess
import os

#os.mkdir("./gff3")

def run_prodigal_phage(inputfasta, out_gene, out_prot):
    to_run="prodigal -i "+inputfasta+" -o "+out_gene+" -a "+out_prot+" -p meta"
    subprocess.call(to_run.split(" "))
    
def run_blastp(inputfasta="input.fasta", output_file="blast.out", database="databast.fasta", evalue="0.001"):
    to_run="blastp -db "+database+" -query "+inputfasta+"-evalue "+evalue+" -outfmt 6 -out "+output_file
    subprocess.call(to_run.split(" "))

def run_formatdb(fastafile, protein="yes"):
    dbtype="prot"
    if protein=="no":
        dbtype="nucl"
    to_run="makeblastdb -in "+fastafile+" -dbtype "+dbtype
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
    
def find_best_hit(gene_id, dict_list):
    evals=1
    annotation=""
    best_hit=""
    for i in range(0, len(first_look)):
        if gene_id in dict_list[i].keys():
            hit=dict_list[i][gene_id]
            #print hit[1]
            if float(hit[1])<evals:
                evals=float(hit[1])
                best_hit=dict_names[i]
                annotation=hit[-1]
            #print(dict_names[i]+"\t"+hit[1]+"\t"+hit[-1])
        
    #print("best annotation for"+gene_id+" is from "+best_hit+" with e-value "+str(evals)+" and annotation of "+annotation)
    return [annotation, best_hit]

#considers hits to more informative databases before less informative databases
#dict_lists are lists of blast_dict tables and dl_names are the names of the dicts in the same order
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
            es.append(float(hit[1]))
            names.append(dl1_names[i])
    if len(hits)>0:
        best_annotation=[hits[es.index(min(es))],names[es.index(min(es))]]
    else:
        for i in range(0, len(dict_list2)):
            if gene_id in dict_list2[i].keys():
                hit=dict_list2[i][gene_id]
                hits.append(hit[-1])
                es.append(hit[1])
                names.append(dl2_names[i])
        if len(hits)>0:
            best_annotation=[hits[es.index(min(es))],names[es.index(min(es))]]
        else:
            best_annotation=["",""]

    #print("best annotation for"+gene_id+" is from "+best_hit+" with e-value "+str(evals)+" and annotation of "+annotation)
    return best_annotation 

#create DB dicts for BLAST file analysis:

#ACLAME:
aclame=open("./databases/DB_Info/aclame/aclame_proteins_all_0.4.tab").readlines()
aclame_dict={}

for line in aclame[2:-1]:
    protein=line.split("\t")[0]
    ncbi_ann=line.split("\t")[2]
    aclame_dict[protein]=ncbi_ann

##COG needs two dbs:
cogs=open("./databases/DB_Info/COG/cog2003-2014.csv").readlines()         #all COG sequences and COG groups
cogs2=open("./databases/DB_Info/COG/cognames2003-2014.tab").readlines()   #COG group definitions/functions

cog_dict={}
cog_defs={}

#dict from gi to COG:
for line in cogs:
    gi=line.split(",")[0]
    cog=line.split(",")[6]
    cog_dict[gi]=cog
#dict from COG to function definition:
for line in cogs2:
    cog=line.split("\t")[0]
    func=line.split("\t")[2]
    cog_defs[cog]=func

##Pfam needs two DBs as well:
pfams=open("./databases/DB_Info/PFam/Pfam-A.titles.txt").readlines()   #ID all pfam sequences in BLAST db
pfams2=open("./databases/DB_Info/PFam/Pfam-A.clans.tsv").readlines()   #Matches IDs to "clans" with functions
pfam_dict={}
pfam_defs={}

#sequence to pfam:
for line in pfams:
    seq=line.split(" ")[0].replace(">","")
    pfam=line.split(" ")[2].split(";")[0]
    pfam_dict[seq]=pfam
#pfam to function:
for line in pfams2:
    pfam=line.split("\t")[0]
    function=line.split("\t")[4]
    pfam_defs[pfam]=function
    
##CAMERA Viral Proteins annotations are complicated; extracting definitions from sequence titles:
cvp=open("./databases/DB_Info/CVP/CVP_titles.txt").readlines()

cvp_dict={}

headers=[]

for line in cvp:
    defs={}
    line=line.replace(">","")
    annotation=""
    info=line.split("/")[1:]
    ID=line.split("/")[0].replace(" ","")
    for i in info: 
        if "=" in i:
            defs[i.split("=")[0]]=i.split("=")[1]
    if "DESCRIPTION" in defs.keys():
        annotation=defs["DESCRIPTION"]
    elif "definition" in defs.keys():
        annotation=defs["definition"]
    
    cvp_dict[ID]=annotation


    
#modules for blast file processing
def get_digits(prodfile):
    prod=open(prodfile).readlines()
    digits=len(str(len(prod)/2))
    return digits


def get_locus_tag(line, digits, phage):
    query=line.split("\t")[0].split(" ")[0]
    number=query.split("_")[-1]
    z="0"*(digits-len(number))
    return "NVP"+phage.replace(".","")+"_"+z+number

#create a dict of BLAST results for each database BLAST:

#ACLAME:
def create_aclame_blast_dict(blast, aclame_dict=aclame_dict, digits=4):
    blast=open(blast).readlines()
    aclame_blast_dict={}

    for line in blast:
        locus_tag=get_locus_tag(line, digits, phage)
        hit=line.split("\t")[1]
        e=line.split("\t")[-2]
        bitscore=line.split("\t")[-1].replace("\n","")
        annotation=aclame_dict[hit]
        info=[hit, e, bitscore, annotation]
        if locus_tag not in aclame_blast_dict.keys():
            aclame_blast_dict[locus_tag]=info
    return aclame_blast_dict
    #print(locus_tag+"\t"+hit+"\t"+e+"\t"+bitscore+"\t"+aclame_dict[hit])

#COGs:
def create_cog_blast_dict(blast, cog_dict=cog_dict, cog_defs=cog_defs, digits=4):
    blast=open(blast).readlines()
    cog_blast_dict={}

    for line in blast:
        locus_tag=get_locus_tag(line, digits, phage)
        hit=line.split("\t")[1]
        e=line.split("\t")[-2]
        bitscore=line.split("\t")[-1].replace("\n","")
        cog=cog_dict[(hit.split("|")[1])]
        func=cog_defs[cog].replace("\n","")
        info=[hit, e, bitscore, cog, func]
        if locus_tag not in cog_blast_dict.keys():
            cog_blast_dict[locus_tag]=info
    return cog_blast_dict

#PFam:
def create_pfam_blast_dict(blast, pfam_dict=pfam_dict, pfam_defs=pfam_defs, digits=4):
    blast=open(blast).readlines()
    pfam_blast_dict={}
    
    for line in blast:
        locus_tag=get_locus_tag(line, digits, phage)
        hit=line.split("\t")[1]
        e=line.split("\t")[-2]
        bitscore=line.split("\t")[-1].replace("\n","")
        pfam=pfam_dict[hit].split(".")[0]
        function=pfam_defs[pfam].replace("\n","")
        info=[hit, e, bitscore, pfam, function]
        if locus_tag not in pfam_blast_dict.keys():
            pfam_blast_dict[locus_tag]=info
    return pfam_blast_dict

##CVP:

def create_cvp_blast_dict(blast, cvp_dict=cvp_dict, digits=4):
    blast=open(blast).readlines()
    cvp_blast_dict={}

    for line in blast:
        locus_tag=get_locus_tag(line, digits, phage)
        hit=line.split("\t")[1]
        e=line.split("\t")[-2]
        bitscore=line.split("\t")[-1].replace("\n","")
        cvp=hit.replace(" ","")
        func=cvp_dict.get(cvp,"")
        info=[hit, e, bitscore, cvp, func]
        if locus_tag not in cvp_blast_dict.keys():
            cvp_blast_dict[locus_tag]=info
    return cvp_blast_dict


def write_cds_gff3(gff, phage, out):
    out=open(out, "w")
    prod=open(gff).readlines()
    digits=len(str(len(prod)/2))  ##for assigning a gene number with appropriate number of zeros preceding it
    
    OGs={"Pfam":pfam_blast_dict, "COG":cog_blast_dict}
    first_looks=[pfam_blast_dict, aclame_blast_dict, cog_blast_dict]
    flnames=["pfam","aclame","cog"]
    second_look=[cvp_blast_dict]
    slnames=["cvp"]
    
    SeqID=prod[0].split(";")[2].split("=")[1].replace('"','')
    
    for i in range(2,len(prod)-1,2):
        loc=prod[i]

        info=prod[i+1]

        if len(loc.split())==2:
            if "complement" in prod[i].split()[1]:
                strand="-"
                start=loc.split("..")[1].replace(")\n","")
                stop=loc.split("(")[1].split("..")[0]
            else:
                strand="+"
                start=loc.split()[1].split("..")[0]
                stop=loc.split()[1].split("..")[1].replace("\n","")
        number=info.split(";")[0].split("=")[2].split("_")[1]
        z="0"*(digits-len(number))
        t="NVP"+phage.replace(".","")+"_"+z+number
        col9="ID="+t
        best_hits=find_best_hit2(t, dict_list1=first_looks, dl1_names=flnames, dict_list2=second_look, dl2_names=slnames)
        Name=best_hits[0]
        if len(Name)==0:
            col9+=', Name=hypothetical protein'
        else:
            col9+=", Name="+Name.replace('"','')
            col9+=", note=annotation from "+best_hits[1]

        for d in OGs.keys():
            db=OGs[d]
            if t in db.keys():
                col9+=', Ontology_term="'+d+":"+db[t][-2]+'"'
        col9+="\n"
        out.write(SeqID+"\tprod\tCDS\t"+start+"\t"+stop+"\t.\t"+strand+"\t0\t"+col9)
    out.close()

import glob

phages=glob.glob("./genomes/*")
phage_list=[]

for p in phages:
    phage_list.append(p.replace("final.fasta","").replace("./genomes/",""))
                      
for phage in phage_list:
    gff="./genes/"+phage+"gene"
    digits=get_digits(gff)
    aclame_blast_dict=create_aclame_blast_dict(blast="./blasts/aclame/"+phage+"vs.aclame.out", digits=digits)
    cog_blast_dict=create_cog_blast_dict(blast="./blasts/cogs_2003-2014/"+phage+"vs.cogs_2003-2014.out", digits=digits)
    pfam_blast_dict=create_pfam_blast_dict(blast="./blasts/Pfam/"+phage+"vs.Pfam.out", digits=digits)
    cvp_blast_dict=create_cvp_blast_dict(blast="./blasts/CVP/"+phage+"vs.CVP.out", digits=digits)
    write_cds_gff3(gff, phage, "./gff3/"+phage+"cds.gff3")
