#!usr/bin/python

import sqlite3
import cPickle as pickle

aclame_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/aclame_dict.p","rb"))
cog_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/cog_dict.p","rb"))
cog_defs=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/cog_def.p","rb"))
pfam_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/pfam_dict.p","rb"))
pfam_defs=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/pfam_def.p","rb"))
cvp_dict=pickle.load(open("/nobackup1/jbrown/annotation/databases/pickled_dicts/cvp_dict.p","rb"))

dd={"aclame":aclame_dict, "cog1":cog_dict, "cog2":cog_defs, "pfam1":pfam_dict,"pfam2":pfam_defs,"cvp":cvp_dict}

conn=sqlite3.connect('/pool001/jbrown/blast_db.sqlite')
c=conn.cursor()

for d in dd.keys():
    db=dd[d]
    if "2" in d:
        c.execute("CREATE TABLE "+d+" (OG primary key, function)")
        for i in db.keys():
            OG=i
            function=db[i].replace("'","''")
            c.execute("INSERT INTO "+d+" VALUES ('"+OG+"','"+function+"')")
    elif "1" in d:
        c.execute("CREATE TABLE "+d+" (ID primary key, OG)")
        for i in db.keys():
            ID=i
            OG=db[i]
            c.execute("INSERT INTO "+d+" VALUES ('"+ID+"','"+OG+"')")
    else: 
        c.execute("CREATE TABLE "+d+" (ID primary key, function)")
        for i in db.keys():
            ID=i
            function=db[i].replace("'","''")
            c.execute("INSERT INTO "+d+" VALUES ('"+ID+"','"+function+"')")

conn.commit()
conn.close()