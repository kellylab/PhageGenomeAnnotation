#!/usr/bin/python

from __future__ import print_function
import subprocess
import glob
import os
import os.path as op
import sys
import click

'''
run_crt.py phage-list --outdir /nobackup1/jbrown/newmu/crt/ --genomedir /nobackup1/jbrown/newmu/genomes 1.028.O
run_crt.py from-genomedir --outdir /nobackup1/jbrown/newmu/crt/ --genomedir /nobackup1/jbrown/newmu/genomes
'''

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    '''run crt for annotation'''
    pass


def run_crt(input_fasta, output, path_to_crt="/home/jbrown/programs/CRT1.2-CLI.jar"):
    args="java -cp "+path_to_crt+" crt -minNR 2 "+input_fasta+" "+output
    subprocess.call(args.split(" "))
    return output


def prep_outdir(out_path):
    if op.exists(out_path) == False:
        os.mkdir(out_path)


def run_crts(phage_names, outdir, genome_dir, path_to_crt="/home/jbrown/programs/CRT1.2-CLI.jar"):
    prep_outdir(outdir)
    for p in phage_names:
        fasta = glob.glob(op.join(genome_dir, "{}*.f*a".format(p)))[0]
        print("found the fasta file: {}".format(fasta))
        output=op.join(outdir, '{}.crt'.format(p))
        output = run_crt(fasta, output, path_to_crt)
        print("CRT finished for {}.  Written to {}".format(p, output))
    return outdir


@cli.command("phage-list", short_help="provide space separated list of phages to blast")
@click.argument('phage-list', nargs=-1)
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--genomedir', help="where phage genomes are found")
@click.option('--path-to-crt', help="where crt executable is found", default='/home/jbrown/programs/CRT1.2-CLI.jar')
def crt_from_list(phage_list, outdir, genomedir, path_to_crt):
    run_crts(phage_list, outdir, genomedir, path_to_crt)


@cli.command("from-genomedir", short_help="provide space separated list of phages to blast")
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--genomedir', help="where phage genomes are found")
@click.option('--path-to-crt', help="where crt executable is found", default='/home/jbrown/programs/CRT1.2-CLI.jar')
def crt_from_dir(outdir, genomedir, path_to_crt):
    fastas = glob.glob(op.join(genomedir, "*.f*a"))
    phage_list = [".".join(op.basename(i).split(".")[:3]) for i in fastas]
    run_crts(phage_list, outdir, genomedir, path_to_crt)


if __name__=='__main__':
    cli()
