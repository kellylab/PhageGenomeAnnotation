#!/usr/bin/python
from __future__ import print_function
import subprocess
import glob
import os
import os.path as op
import click

'''
phage_ublast_all.py phage-list --outdir /nobackup1/jbrown/newmu/blasts/ --proteindir /nobackup1/jbrown/newmu/proteins/ 1.028.O
'''

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    '''run ublasts for annotation'''
    pass

udbs=["aclame.udb", "cogs_2003-2014.udb", "CVP.udb", "eggnog4.udb", "Pfam.udb", "pog.udb", "tara.translated.udb","kegg.reduced.fasta"]
dbnames=["aclame","cogs_2003-2004","CVP","eggnog","Pfam","pog","tara.translated","kegg"]

def run_ublastp(fastafile, out, udb, evalue, ublast_path='/home/sbiller/usearch7.0.1090_i86linux64'):
    to_run = ublast_path + " -ublast "+fastafile+" -db "+udb+" -evalue "+evalue+" -accel 0.5 -strand plus -blast6out "+out
    subprocess.call(to_run.split(" "))

def prep_outdirs(outdir = "/nobackup1/jbrown/annotation/blasts/", dbnames=dbnames):
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    for dbname in dbnames:
        if os.path.exists(op.join(outdir,dbname)) == False:
            os.mkdir(op.join(outdir,dbname))

def run_ublasts(phage_list, outdir, databasedir, proteindir, ublast_path='/home/sbiller/usearch7.0.1090_i86linux64', evalue='1e-5'):
    prep_outdirs(outdir)
    for p in phage_list:
        phage_proteins = glob.glob(op.join(proteindir, "{i}*faa".format(i=p)))[0]
        for i in range(0, len(udbs)):
            udb = op.join(databasedir, udbs[i])
            dbname = dbnames[i]
            print("running a ublast comparison of {} against {}".format(phage_proteins, udb))
            out_blast = op.join(outdir, dbname, p + "vs."+dbname+".out")
            if op.exists(out_blast):
                continue
            else:
                run_ublastp(fastafile=phage_proteins, udb=udb, out=op.join(outdir, dbname, p + "vs."+dbname+".out"), evalue=evalue, ublast_path=ublast_path)
    return outdir


@cli.command("phage-list", short_help="provide space separated list of phages to blast")
@click.argument('phage_list', nargs=-1)
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--databasedir', help="where to find the blast databases", default='/nobackup1/jbrown/annotation/databases')
@click.option('--proteindir', help="where phage proteins are found")
@click.option('--ublast_path', help="where ublast executable found", default='/home/sbiller/usearch7.0.1090_i86linux64')
def main(phage_list, outdir, databasedir, proteindir):
    run_ublasts(phage_list, outdir, databasedir, proteindir)


@cli.command("from-proteindir", short_help="find proteomes within indicated protein directory")
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--databasedir', help="where to find the blast databases", default='/nobackup1/jbrown/annotation/databases')
@click.option('--proteindir', help="where phage proteins are found")
@click.option('--nvp', help='indicate whether protein files follow nvp naming scheme', default=True, show_default=True)
def main_alldir(outdir, databasedir, proteindir, nvp):
    fastas = glob.glob(op.join(proteindir, "*.f*a"))
    if nvp is True:
        phage_list = [".".join(op.basename(i).split(".")[:3]) for i in fastas]
    else:
        phage_list = [i.split(".")[0] for i in fastas]

    run_ublasts(phage_list, outdir, databasedir, proteindir)


if __name__ == '__main__':
    cli()
