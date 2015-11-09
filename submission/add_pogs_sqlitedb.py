#!usr/bin/python
import sqlite3

conn=sqlite3.connect('/pool001/jbrown/blast_db.sqlite')
c=conn.cursor()

with open("/nobackup1/jbrown/annotation/databases/pogs.txt") as infile: 
    c.execute("CREATE TABLE "+"pog"+" (ID, OG, function, phy1, phy2)")
    for li in infile:
        ID= li.split(":")[3]
        pog=li.split(":")[0]
        function=li.split(":")[5].split("|")[-1].replace("'","''")
        try:
            phy1=li.split(":")[6].replace("'","''")
        except:
            phy1="none found"
        try:
            phy2=li.split(":")[7].replace("'","''")
        except:
            phy2="none found"
        
        c.execute("INSERT INTO pog VALUES ('%s','%s','%s','%s','%s')" %(ID, pog, function, phy1,phy2))
conn.commit()
conn.close()