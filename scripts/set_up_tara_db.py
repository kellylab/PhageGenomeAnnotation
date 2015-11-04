#!usr/bin/python

import sqlite3

conn=sqlite3.connect('tara_db.sqlite')
c=conn.cursor()
c.execute('''CREATE TABLE taratbl
            (ID primary key, gene, egg, ko, kfunc)''')
with open("./OM-RGC_seq.release.tsv") as infile:
    for r in infile:
        vec=r.replace('"','').split("\t")
        ID=vec[0]
        gene=vec[1]
        egg=vec[2]
        ko=vec[3]
        kfunc=vec[4]
        c.execute("INSERT INTO taratbl VALUES ('"+ID+"','"+gene+"','"+egg+"','"+ko+"','"+kfunc+"')")
conn.commit()
conn.close()