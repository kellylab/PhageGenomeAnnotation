#!/usr/bin/python
import glob

gff3s=glob.glob("./gff3/*")
#gff3s=["1.161.O.cds.gff3"]

out=open("cds.annotation.info1.txt","w")
out.write("phage\ttotal_orfs\ttotal_cogs\ttotal_pfam\thypotheticals\n")

for g in gff3s:
    phage=g.split("/")[-1].replace("cds.gff3","")
    cog_ct=0
    pfam_ct=0
    unknowns=0
    ann=open(g).readlines()
    total=len(ann)
    for line in ann:
        if "COG" in line:
            cog_ct+=1
        if "Pfam" in line:
            pfam_ct+=1
        if "hypothetical" in line:
            unknowns+=1
    out.write(phage+"\t"+str(total)+"\t"+str(cog_ct)+"\t"+str(pfam_ct)+"\t"+str(unknowns)+"\n")
out.close()
   
    
