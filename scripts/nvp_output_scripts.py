#!usr/bin/python
from nvp_blast_processing import *

import os
from pyfaidx import Fasta

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
            stop=loc.split("..")[1].replace(")\n","")
            start=loc.split("(")[1].split("..")[0]
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

gene_id="NVP1161O_169"

def id_blast_hits(gene_id, blast_dict, ann1):
    hits=[]
    es=[]
    names=[]
    for i in range(0, len(ann1)):
        if gene_id in blast_dict[ann1[i]].keys():
            hit=blast_dict[ann1[i]][gene_id]
            hits.append(hit[-1])
            es.append(float(hit[3]))
            names.append(ann1[i])
    return [hits, es, names]


def find_best_hit2(gene_id, blast_dict, ann1=["cog","pfam","aclame","kegg","egg","pog"]):
    evals=1
    annotation=""
    best_hit=""
    ids=id_blast_hits(gene_id, blast_dict, ann1)
    hits=ids[0]
    es=ids[1]
    names=ids[2]
    if len(ids[0])>0:
        best_annotation=[hits[es.index(min(es))],names[es.index(min(es))]]
        if re.search(r"none|^NA|[\d]{1,5}\.\.[\d]{1,5}|hypothetical", best_annotation[0]):
            ann1.remove(best_annotation[1])
            find_best_hit2(gene_id, blast_dict, ann1=ann1)
        else:
            return best_annotation
    else:
        best_annotation=["",""]
    return best_annotation 

def gff3_header(prod):
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    return Sequence+"\n"

#merge BLAST results into one gff3
    
def cds_blast_annotations_to_gff3(phage, cov_thresh=75):
    #prodigal and fasta files:
    prod=prod_path+phage+"gene"
    faa=faa_path+phage+"faa"
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    
    ogs=["kegg","pfam","cog","egg", "pog"]
    annotes=["aclame","tara","cvp"]
    
    blast_dict=load_blast_files(phage, cov_thresh=cov_thresh)
    
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
        best_hits=find_best_hit2(locus_tag, blast_dict)

        #establish name:
        Name=best_hits[0]
        if len(Name)==0:
            col9+='; Name=hypothetical protein'
        else:
            col9+="; Name="+Name.replace('"','')

        #Add OG annotations:
        for d in range(0, len(ogs)):
            og_dict=blast_dict[ogs[d]]
            if locus_tag in og_dict.keys():
                col9+='; Ontology_term="'+ogs[d]+":"+og_dict[locus_tag][-2]+'"'

        #Add db closest hits to notes
        for d in range(0, len(annotes)):
            annote_dict=blast_dict[annotes[d]]
            if locus_tag in annote_dict.keys():
                col9+='; note="'+annotes[d]+"_best_match:"+annote_dict[locus_tag][-2]+'"'
        out+=Sequence+"\t"+"prod"+"\t"+"CDS"+"\t"+start+"\t"+stop+"\t"+"."+"\t"+strand+"\t"+"0"+"\t"+col9+"\n"
    return out

def CRISPR_gff3(phage):
    crt_output="/nobackup1/jbrown/annotation/crt/"+phage+"crt"
    crtout=open(crt_output).readlines()
    name=phage
    SeqID=crtout[0].split()[1]
    out=""
    for line in crtout:
        if line.startswith("CRISPR"):
            vec=line.split()
            number=vec[1]
            start=vec[3]
            stop=vec[5]
            ID="NVP"+name.replace(".","")+"_CRISPR-like_"+number
            out+=SeqID+"\t"+"crt"+"\tputative CRISPR feature\t%s\t%s\t.\t.\t.\tID=%s" % (start, stop, ID)
            out+=", note=CRISPR region\n"
    return out

def tRNA_scan_to_gff3(phage):
    tRNAScanSE_file="/nobackup1/jbrown/annotation/trna/%strna" %phage
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
            anticodon=l[5]
            SeqID=l[0]
            col9="ID="+locus_tag+", aa="+aa+", anticodon="+anticodon
            out=SeqID+"\t"+"tRNAScanSE"+"\t"+"tRNA"+"\t"+start+"\t"+stop+"\t"+l[-1].replace("\n","")+"\t"+strand+"\t"+"0"+"\t"+col9+"\n"
            tanns+=out
        return tanns
    else:
        print "no tRNAs found in genome"
        return ""
    
#put them all together:
def write_gff3_file(phage, output_file, cov_thresh=75):
    prod=prod_path+phage+".gene"
    faa=faa_path+phage+".faa"
    genomic_fasta="/nobackup1/jbrown/annotation/genomes/%sfinal.fasta" % phage
    
    out=open(output_file,"w")
    #out.write(gff3_header(phage+"gene"))
    out.write(cds_blast_annotations_to_gff3(phage, cov_thresh=cov_thresh))
    
    
    out.write(tRNA_scan_to_gff3(phage))
    
    crt_output=phage+"crt"
    out.write(CRISPR_gff3(phage))
    
    out.close()

    
#function to query sqlite tara db... not going to use because tara project used very weak
#evalues to assign function/OGs to sequences within their library
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

def cds_blast_annotations_to_table(phage, cov_thresh=75):
    #prodigal and fasta files:
    prod=prod_path+phage+"gene"
    faa=faa_path+phage+"faa"
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    
    ogs=["kegg","pfam","cog","egg"]
    annotes=["aclame","tara","cvp"]
    
    blast_dict=load_blast_files(phage, cov_thresh=cov_thresh)
    
    out=""  #set up string to write to
    #write gff3 lines from prodigal files and blast dicts:
    out=""
    #out+="genome\tlocus_tag\ttype\tstart\tstop\tstrand\tbest_hit_annotation\t"
    #out+="kegg\tpfam\tcog\tegg\taclame\ttara\tcvp\n"
    
    #run through annotations of each prodigal-identified CDS:
    prod=open(prod).readlines()
    digits=get_digits(faa)
    
    #write lines from prodigal files and blast dicts:
    for i in range(2,len(prod)-1,2):
        out+=Sequence+"\t"
        
        
        coords=get_prod_cds_info(i, prod, digits, phage)
        locus_tag=coords[0]
        start=coords[1]
        stop=coords[2]
        strand=coords[3]
        
        out+=locus_tag+"\tcds\t"+start+"\t"+stop+"\t"+strand+"\t"
        
        #ID best hit:
        best_hits=find_best_hit2(locus_tag, blast_dict)

        #establish name:
        Annotation=best_hits[0]
        if len(Annotation)==0:
            Annotation='hypothetical protein'
        else:
            Annotation=Annotation.replace('"','')
        out+=Annotation+"\t"
        
        #Add OG annotations:
        for d in range(0, len(ogs)):
            og_dict=blast_dict[ogs[d]]
            if locus_tag in og_dict.keys():
                out+=og_dict[locus_tag][-2]+"\t"
            else:
                out+="NA\t"

        #Add db closest hits to notes
        for d in range(0, len(annotes)):
            annote_dict=blast_dict[annotes[d]]
            if locus_tag in annote_dict.keys():
                out+=annote_dict[locus_tag][-2]+'\t'
            else:
                out+="NA\t"
        out+="\n"

    return out

def kegg_egg_pog_tbl(phage, cov_thresh=75):
    #prodigal and fasta files:
    prod=prod_path+phage+"gene"
    faa=faa_path+phage+"faa"
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')
    
    ogs=["kegg","egg","pog"]
    og_lens={"kegg":7,"egg":7,"pog":9}
    blast_dict=load_kegg_egg_pog_blast(phage, cov_thresh=cov_thresh)
    
    out=""  #set up string to write to

    
    #run through annotations of each prodigal-identified CDS:
    prod=open(prod).readlines()
    digits=get_digits(faa)
    
    #write lines from prodigal files and blast dicts:
    for i in range(2,len(prod)-1,2):
        #col1
        out+=Sequence+"\t"
        
        
        locus_tag=get_prod_cds_info(i, prod, digits, phage)[0]
        #col2
        out+=locus_tag+"\t"
        
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