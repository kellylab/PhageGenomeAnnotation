#!usr/bin/python
import sqlite3

'''
#skipping because already in db:
conn=sqlite3.connect('/pool001/jbrown/blast_db.sqlite')
c=conn.cursor()

with open("/nobackup1/jbrown/annotation/databases/DB_Info/eggNOG/NOG.members.tsv") as infile: 
    c.execute("CREATE TABLE "+"egg1"+" (ID, OG, category)")
    for li in infile:
        IDs= li.split("\t")[-1]
        cog=li.split("\t")[1]
        func_cat=li.split("\t")[4]
        for i in IDs.split(","):
            c.execute("INSERT INTO "+"egg1"+" VALUES ('"+i+"','"+cog+"','"+func_cat+"')")
conn.commit()
conn.close()
'''

conn=sqlite3.connect('/pool001/jbrown/blast_db.sqlite')
c=conn.cursor()

with open("/nobackup1/jbrown/annotation/databases/DB_Info/eggNOG/NOG.annotations.tsv") as infile:
    c.execute("CREATE TABLE "+"egg2"+" (OG primary key, function, category)")
    for l in infile:
        cog=l.split("\t")[1]
        function=l.split("\t")[-1].replace("\n","").replace("'","''")
        func_cat=l.split("\t")[-2]
        c.execute("INSERT INTO egg2 VALUES('%s','%s','%s')" %(cog, function, func_cat))
conn.commit()
conn.close()