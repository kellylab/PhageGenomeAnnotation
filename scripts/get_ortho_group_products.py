 #!usr/bin/python

from Bio import SeqIO
import sqlite3
from collections import Counter
import sys

def set_org_name_dict(orthomcl_codenames):
    'dict of orthomcl taxa_ids and original organism names'
    orgs=open(orthomcl_codenames).readlines()
    org_dict={}
    for o in orgs:
        vec=o.split("\t")
        organism=vec[0].split("/")[-1].replace("fasta","").replace("faa","")
        code=vec[1].replace("\n","")
        org_dict[code]=organism
    return org_dict

class Protein():
    ortho_name_dict=set_org_name_dict( "/nobackup1/jbrown/annotation/orthomcl/host_v_host/host_taxon_codes.txt")
    
    def __init__(self, name):
        self.name=name
    
    def faa_name(self):
        return self.name.split("|")[1]
    
    def locus_tag(self):
        contig=self.faa_name().split("_")[2]
        number=self.faa_name().split("_")[3]
        return "ORF_%s_%s" % (contig, number)
    
    def gbk_file(self):
        code=self.name.split("|")[0]
        long_name=name_dict[code].replace("_cds_prod.","_contigs*_prod.gbk")
        path="/nobackup1/jbrown/vibrio_genomes/gbk/"
        return path+long_name
    
    def product(self):
        conn=sqlite3.connect('/pool001/jbrown/HRX_Vibrio_lt_db.sqlite')
        c=conn.cursor()
        c.execute("SELECT product from products where locus_tag='%s'" % self.faa_name())
        output=c.fetchone()
        if output==None:
            result="not in db"
        else:
            result=output[0]
        conn.close()
        return result
    
def create_group_dict(group_file):
    group_dict={}
    with open(group_file) as infile:
        for line in infile:
            cluster=line.split(":")[0]
            proteins=line.split(":")[1].replace("\n","").split(" ")[1:]
            group_dict[cluster]=proteins
    return group_dict

def get_all_group_annotations(group, group_dict):
    proteins=group_dict[group]
    products=[]
    for p in proteins:
        pobj=Protein(p)
        products.append(pobj.product())
    return products

def get_most_common_group_annotation(group, group_dict):
    proteins=group_dict[group]
    products=[]
    for p in proteins:
        pobj=Protein(p)
        products.append(pobj.product())
    prod_count=Counter(products)
    mc=prod_count.most_common(1)
    return list(mc)[0][0]

def write_best_hits(ortho_group_file, cluster_list, output_file):
    clusters=open(cluster_list).readlines()
    group_dict=create_group_dict(ortho_group_file)
    out=open(output_file,"w")
    out.write("cluster\tmost_common_annotation\tgroup_size\n")
    for c in clusters:
        clust=c.replace("\n","")
        best=get_most_common_group_annotation(clust, group_dict)
        num_prots=len(group_dict[clust])
        out.write("%s\t%s\t%s\n" % (clust, best, num_prots))
    out.close()

cluster_list="/nobackup1/jbrown/annotation/orthomcl/host_v_host/hi_mi_hostprot_phage_pairs.txt"
output="/nobackup1/jbrown/annotation/orthomcl/host_v_host/himi_cluster_annotations_from_gbks.txt"
ortho_group_file="/nobackup1/jbrown/annotation/orthomcl/host_v_host/groups.txt"

write_best_hits(ortho_group_file, cluster_list, output)
              
