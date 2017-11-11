from __future__ import print_function
import subprocess
import glob
import os
import os.path as op
import sys
import click

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    '''run prodigal for annotation'''
    pass

def run_prodigal_phage(inputfasta, out_gene, out_fna, out_prot):
    to_run="prodigal -i "+inputfasta+" -o "+out_gene+" -a "+out_prot+" -d "+out_fna+" -p meta"
    subprocess.call(to_run.split(" "))


def prep_outdir(out_path):
    if op.exists(out_path) == False:
        os.mkdir(out_path)

def run_prodigals(phage_names, outdir, genome_dir):
    prep_outdir(outdir)
    outdirs = [op.join(outdir, i) for i in ['genes', 'proteins','fna']]
    blank = [prep_outdir(d) for d in outdirs]
    for p in phage_names:
        fasta = glob.glob(op.join(genome_dir, "{}*.f*a".format(p)))[0]
        print("found the fasta file: {}".format(fasta))
        out_gene = op.join(outdirs[0],"{}.gene".format(p))
        out_prot = op.join(outdirs[1], "{}.faa".format(p))
        out_fna = op.join(outdirs[2], "{}.fna".format(p))
        output = run_prodigal_phage(fasta, out_gene, out_fna, out_prot)
        print("prodigal finished for {}.  Written to {}".format(p, output))
    return outdirs


@cli.command("phage-list", short_help="provide space separated list of phages to blast")
@click.argument('phage-list', nargs=-1)
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--genomedir', help="where phage genomes are found")
def prodigal_from_list(phage_list, outdir, genomedir):
    run_prodigals(phage_list, outdir, genomedir)


@cli.command("from-genomedir", short_help="provide space separated list of phages to blast")
@click.option('--outdir', help="where to send blast outputs", default="/nobackup1/jbrown/annotation/blasts/")
@click.option('--genomedir', help="where phage genomes are found")
def trnascan_from_dir(outdir, genomedir):
    fastas = glob.glob(op.join(genomedir, "*.f*a"))
    phage_list = [".".join(op.basename(i).split(".")[:3]) for i in fastas]
    run_prodigals(phage_list, outdir, genomedir)


if __name__=='__main__':
    cli()
