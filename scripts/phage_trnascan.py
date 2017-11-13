#!/usr/bin/python
from __future__ import print_function
import subprocess
import glob
import os
import os.path as op
import click

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    '''run trnascan for annotation'''
    pass

def run_trna_scan(input_file, output):
    if os.path.exists(output):
        return ouput
    else:
        args=["tRNAscan-SE", "-o", output, "-G", "-D", input_file]
        subprocess.call(args)

def prep_outdir(out_path):
    if op.exists(out_path) == False:
        os.mkdir(out_path)

def run_trna_scans(phage_names, outdir, genome_dir):
    prep_outdir(outdir)
    for p in phage_names:
        fasta = glob.glob(op.join(genome_dir, "{}*.f*a".format(p)))[0]
        print("found the fasta file: {}".format(fasta))
        output=op.join(outdir, '{}.trna'.format(p))
        output = run_trna_scan(fasta, output)
        print("tRNA Scan finished for {}.  Written to {}".format(p, output))
    return outdir


@cli.command("phage-list", short_help="provide space separated list of phages to blast")
@click.argument('phage-list', nargs=-1)
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--genomedir', help="where phage genomes are found")
def trnascan_from_list(phage_list, outdir, genomedir):
    run_trna_scans(phage_list, outdir, genomedir)


@cli.command("from-genomedir", short_help="find genomes within indicated genome directory")
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--genomedir', help="where phage genomes are found")
@click.option('--nvp', help='indicate whether protein files follow nvp naming scheme', default=True, show_default=True)
def trnascan_from_dir(outdir, genomedir, nvp):
    fastas = glob.glob(op.join(genomedir, "*.f*a"))
    if nvp is True:
        phage_list = [".".join(op.basename(i).split(".")[:3]) for i in fastas]
    else:
        phage_list = [i.split(".")[0] for i in fastas]
    run_trna_scans(phage_list, outdir, genomedir)


if __name__=='__main__':
    cli()
