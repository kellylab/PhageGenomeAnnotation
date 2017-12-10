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

def get_prod_cds_info(i, prod, digits, phage, genome_len):
    loc=prod[i]
    if len(loc.split())==2:
        start_prefix = ""
        stop_prefix = ""
        if "complement" in prod[i].split()[1]:
            strand="-"
            stop=loc.split("..")[1].replace(")\n","")
            start=loc.split("(")[1].split("..")[0]
            if ">" in start:
                start = start.replace(">","")
                start_prefix = ">"
            elif "<" in start:
                start = start.replace("<","")
                start_prefix=">"

            if ">" in stop:
                stop = stop.replace(">","")
                stop_prefix = "<"
            elif "<" in stop:
                stop = stop.replace(">","")
                stop_prefix = "<"
        else:
            strand="+"
            start=loc.split()[1].split("..")[0]
            stop=loc.split()[1].split("..")[1].replace("\n","")
            if "<" in start:
                start = start.replace("<","")
                start_prefix = "<"
            elif ">" in start:
                start = start.replace(">","")
                start_prefix = "<"

            if "<" in stop:
                stop = stop.replace("<","")
                stop_prefix = ">"

            elif ">" in stop:
                stop = stop.replace(">","")
                stop_prefix = ">"

        start=int(start)
        stop=int(stop)
        real_start = int(start)
        real_stop = int(stop)

    if len(start_prefix) == 1:
        if start < 4:
            print("start is less than 4")
            start = 1

        elif start > (genome_len - 4):
            print("start is greater than the genome length")
            print(start)
            print(genome_len - 4)
            start = genome_len

        else:
            print("EXCEPTION FOUND")
    if len(stop_prefix) == 1:
        if stop < 4:
            print("stop is less than 4")
            stop = 1
        elif stop > (genome_len - 4):
            print("stop is greater than the genome length")
            stop = genome_len
        else:
            print('EXCEPTION FOUND')
    start_all = "{}{}".format(start_prefix, start)
    stop_all = "{}{}".format(stop_prefix, stop)
    info=prod[i+1]
    number=info.split(";")[0].split("=")[2].split("_")[1]
    z="0"*(digits-len(number))
    t="NVP"+phage.replace(".","")+"_"+z+number
    return [t, start_all, stop_all, strand, real_start, real_stop]

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

def cds_blast_annotations_to_gff3(phage, prod_path, faa_path, blast_path, cov_thresh=75):
    #prodigal and fasta files:
    prod=glob.glob(op.join(prod_path,phage + "*.gen*"))[0]
    faa=glob.glob(op.join(faa_path, phage+"*.f*a"))[0]
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')

    ogs=["kegg","pfam","cog","egg", "pog"]
    annotes=["aclame","tara","cvp"]

    blast_dict=load_blast_files(phage, prod_path, faa_path, blast_path, cov_thresh=cov_thresh)

    out=""  #set up string to write to

    #run through annotations of each prodigal-identified CDS:
    with open(prod) as ih:
        prod= ih.readlines()
        digits=get_digits(faa)
        seq_len = int(prod[0].split(";")[1].split("=")[1])
        #write gff3 lines from prodigal files and blast dicts:
        for i in range(2,len(prod)-1,2):

            coords=get_prod_cds_info(i, prod, digits, phage, seq_len)
            locus_tag=coords[0]
            start=coords[1]
            stop=coords[2]
            strand=coords[3]

            #set up col9
            col9="ID="+locus_tag

            ### CHANGED:
            if str(start) != str(coords[4]) and strand == "+":
                col9 += "; codon_start=%s" % coords[4]

            if str(stop) != str(coords[5]) and strand == '-':
                new_val = int(stop) - int(coords[5]) + 1
                col9 += "; codon_start=%s" % coords[5]
            ### CHANGED OVER

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

def CRISPR_gff3(phage, crt_path="/nobackup1/jbrown/annotation/crt/"):
    crt_output=op.join(crt_path, phage+".crt")
    with open(crt_output) as cout:
        num = 0
        crtout = cout.readlines()
        name=phage
        SeqID=crtout[0].split()[1]
        out=""
        for i, line in enumerate(crtout):
            if "[" in line:
                num += 1
                vec = line.strip().split("\t")
                start = vec[0]
                repeat_end = int(vec[4].split(",")[0].replace('[ ',""))
                spacer_end = int(vec[4].split(",")[1].replace(']',"").replace(" ",""))
                crispr_start = int(start) + repeat_end
                crispr_stop = int(crispr_start) + spacer_end
                ID="NVP"+name.replace(".","")+"_CRISPR_spacer_"+str(num)
                out+=SeqID+"\t"+"crt"+"\tncRNA\t%s\t%s\t.\t.\t.\tID=%s" % (crispr_start, crispr_stop, ID)
                out+='ncRNA_class="scRNA"; note=CRISPR spacer {num}\n'.format(num=num)
    return out

def tRNA_scan_to_gff3(phage, trna_path="/nobackup1/jbrown/annotation/trna"):
    tRNAScanSE_file=op.join(trna_path,"{}.trna".format(phage))
    if os.path.getsize(tRNAScanSE_file)>0:
        with open(tRNAScanSE_file) as tout:
            t = tout.readlines()
            tanns=""
            for line in t[3:]:
                l=line.split("\t")
                locus_tag="NVP"+l[0].split("_")[1].replace(".","")+"_tRNA_"+l[1]
                start=l[2]
                stop=l[3]
                if start < stop:
                    strand="+"
                else:
                    strand="-"
                if l[4] == "Undet" or l[4] == "Sup":
                    aa = "Xxx"
                else:
                    aa=l[4]
                anticodon=l[5]
                SeqID=l[0]
                col9="ID=tRNA-" + aa + '; note=anticodon:' + anticodon
                out=SeqID+"\t"+"tRNAScanSE"+"\t"+"tRNA"+"\t"+start+"\t"+stop+"\t"+l[-1].replace("\n","")+"\t"+strand+"\t"+"0"+"\t"+col9+"\n"
                tanns+=out
        return tanns
    else:
        print("no tRNAs found in genome for {}".format(phage))
        return ""

#put them all together:
def write_gff3_file(phage, output_file, prod_path, faa_path, genome_path, blast_path, trna_path, crt_path, cov_thresh=75):
    '''
    >>output_file = "/nobackup1/jbrown/newmu/gff3/1.028.O.gff3"
    >>prod_path = "/nobackup1/jbrown/newmu/genes/"
    >>faa_path = "/nobackup1/jbrown/newmu/proteins/"
    >>genome_path = "/nobackup1/jbrown/newmu/genomes/"
    >>blast_path = "/nobackup1/jbrown/newmu/blasts/"
    >>trna_path = "/nobackup1/jbrown/newmu/trna/"
    >>crt_path = "/nobackup1/jbrown/newmu/crt/"
    >> write_gff3_file('1.028.O', output_file, prod_path, faa_path, genome_path, blast_path, trna_path, crt_path, cov_thresh=75)
    '''
    print("prod dir: {}".format(prod_path))
    print("looking for: {}".format(op.join(prod_path,phage + "*.gen*")))
    prod=glob.glob(op.join(prod_path,phage + "*.gen*"))[0]
    faa=glob.glob(op.join(faa_path, phage+"*.f*a"))[0]
    genomic_fasta=op.join(genome_path, "%s.final.fasta" % phage)

    with open(output_file,"w") as out:
    #cds_blast_annotations_to_gff3(phage, prod_path, faa_path, blast_path, cov_thresh=75):
        print("identifying annotations based on blast results for {}".format(phage))
        out.write(cds_blast_annotations_to_gff3(phage, prod_path, faa_path, blast_path, cov_thresh=cov_thresh))
        #tRNA_scan_to_gff3(phage, trna_path="/nobackup1/jbrown/annotation/trna"):
        print("identifying trna annotations for {}".format(phage))
        out.write(tRNA_scan_to_gff3(phage, trna_path))
        out.write(CRISPR_gff3(phage, crt_path))


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

def cds_blast_annotations_to_table(phage, prod_path, faa_path, blast_path, cov_thresh=75):
    #prodigal and fasta files:
    prod=glob.glob(op.join(prod_path,phage + "*.gen*"))[0]
    faa=glob.glob(op.join(faa_path, phage+"*.f*a"))[0]
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')

    ogs=["kegg","pfam","cog","egg"]
    annotes=["aclame","tara","cvp"]
    #load_blast_files(phage, prod_path, faa_path, blast_path, cov_thresh=75):
    blast_dict=load_blast_files(phage, prod_path, faa_path, blast_path, cov_thresh=cov_thresh)

    out=""  #set up string to write to
    #write gff3 lines from prodigal files and blast dicts:
    out=""
    #out+="genome\tlocus_tag\ttype\tstart\tstop\tstrand\tbest_hit_annotation\t"
    #out+="kegg\tpfam\tcog\tegg\taclame\ttara\tcvp\n"

    #run through annotations of each prodigal-identified CDS:
    with open(prod) as ih:
        prod = ih.readlines()
        digits=get_digits(faa)
        seq_len = int(prod[0].split(";")[1].split("=")[1])

        #write lines from prodigal files and blast dicts:
        for i in range(2,len(prod)-1,2):
            out+=Sequence+"\t"
            coords=get_prod_cds_info(i, prod, digits, phage, seq_len)
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
    prod=glob.glob(op.join(prod_path,phage + "*.gen*"))[0]
    faa=glob.glob(op.join(faa_path, phage+"*.f*a"))[0]
    Sequence=open(prod).readlines()[0].split(";")[2].split("=")[1].replace('"','')

    ogs=["kegg","egg","pog"]
    og_lens={"kegg":7,"egg":7,"pog":9}
    blast_dict=load_kegg_egg_pog_blast(phage, cov_thresh=cov_thresh)

    out=""  #set up string to write to


    #run through annotations of each prodigal-identified CDS:
    prod=open(prod).readlines()
    digits=get_digits(faa)
    seq_len = int(prod[0].split(";")[1].split("=")[1])

    #write lines from prodigal files and blast dicts:
    for i in range(2,len(prod)-1,2):
        #col1
        out+=Sequence+"\t"


        locus_tag=get_prod_cds_info(i, prod, digits, phage, seq_len)[0]
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
