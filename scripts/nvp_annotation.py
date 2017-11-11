from __future__ import print_function
import click
import os.path as op
import os
import glob

from phage_ublast import run_ublasts
from phage_prodigal import run_prodigals
from phage_trnascan import run_trna_scans, prep_outdir
from phage_crt import run_crts
from nvp_output_scripts import write_gff3_file

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    '''functions to output gff3s'''
    pass

def write_gff3s(phage_list, output_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_threshold=75):
    prep_outdir(output_dir)
    for phage in phage_list:
        output_file = op.join(output_dir, '{}.gff3'.format(phage))
        write_gff3_file(phage, output_file, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_thresh=75)
    return outdir


@cli.command("phage-list-runall", short_help="provide space separated list of phages to blast")
@click.argument('phage-list', nargs=-1)
@click.option('--outdir', help="where to send blast outputs")
@click.option('--genome-dir', help="where to find genomic contigs in fasta format")
@click.option('--blast_databasedir', help="where to find the blast databases", default='/nobackup1/jbrown/annotation/databases')
@click.option('--ublast_path', help="where ublast executable found", default='/home/sbiller/usearch7.0.1090_i86linux64')
@click.option('--ublast_evalue', help='evalue to use for ublast', default='1e-5')
@click.option('--path-to-crt', help='location of crt executable', default="/home/jbrown/programs/CRT1.2-CLI.jar")
def run_all(phage_list, genome_dir, outdir, blast_databasedir, ublast_path, ublast_evalue):
    print("Running prodigal now")
    prod_dir, prot_dir, fna_dir = run_prodigals(phage_list, outdir, genome_dir)

    blast_dir = op.join(outdir, "blasts")
    print("Running blasts now")
    blast_dir = run_ublasts(phage_list, blast_dir, blast_databasedir, prot_dir, ublast_path, ublast_evalue)

    trna_dir = op.join(outdir, 'trna')
    print("Running trna scan now")
    trna_dir = run_trna_scans(phage_list, trna_dir, genome_dir)

    crt_dir = op.join(outdir, 'crt')
    print("Looking for CRISPRs using crt")
    crt_dir = run_crts(phage_list, crt_dir, genome_dir, path_to_crt=path_to_crt)

    gff_dir = op.join(outdir, 'gff3')
    print("Writing results to gff3")
    gff_dir = write_gff3s(phage_list, gff_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_threshold=75)


@cli.command("genome-dir-runall", short_help="provide space separated list of phages to blast")
@click.option('--outdir', help="where to send blast outputs")
@click.option('--genome-dir', help="where to find genomic contigs in fasta format")
@click.option('--blast_databasedir', help="where to find the blast databases", default='/nobackup1/jbrown/annotation/databases')
@click.option('--ublast_path', help="where ublast executable found", default='/home/sbiller/usearch7.0.1090_i86linux64')
@click.option('--ublast_evalue', help='evalue to use for ublast', default='1e-5')
@click.option('--path-to-crt', help='location of crt executable', default="/home/jbrown/programs/CRT1.2-CLI.jar")
def run_all(phage_list, genome_dir, outdir, blast_databasedir, ublast_path, ublast_evalue):
    fastas = glob.glob(op.join(genome_dir, "*.f*a"))
    phage_list = [".".join(op.basename(i).split(".")[:3])+"." for i in fastas]

    print("Running prodigal now")
    prod_dir, prot_dir, fna_dir = run_prodigals(phage_list, outdir, genome_dir)

    blast_dir = op.join(outdir, "blasts")
    print("Running blasts now")
    blast_dir = run_ublasts(phage_list, blast_dir, blast_databasedir, prot_dir, ublast_path, ublast_evalue)

    trna_dir = op.join(outdir, 'trna')
    print("Running trna scan now")
    trna_dir = run_trna_scans(phage_list, trna_dir, genome_dir)

    crt_dir = op.join(outdir, 'crt')
    print("Looking for CRISPRs using crt")
    crt_dir = run_crts(phage_list, crt_dir, genome_dir, path_to_crt=path_to_crt)

    gff_dir = op.join(outdir, 'gff3')
    print("Writing results to gff3")
    gff_dir = write_gff3s(phage_list, gff_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_threshold=75)


@cli.command("write-gff3s", short_help='write gff3 files given all appropriate input directories')
@click.option("--genome-dir", help="location of genomic contigs in fasta format")
@click.option("--outdir", help="parent directory of all other results directories, assumes that trnascan, crt and blasts have already been run and have default naming structure")
def write_from_genomedir(genome_dir, outdir):
    blast_dir = op.join(outdir, "blasts")
    trna_dir = op.join(outdir, 'trna')
    crt_dir = op.join(outdir, 'crt')
    prod_dir = op.join(outdir, 'genes')
    prot_dir = op.join(outdir, 'proteins')

    for i in [blast_dir, trna_dir, crt_dir]: assert op.exists(i), "Please make sure that {i} exists with a result per genome".format(i=i)

    fastas = glob.glob(op.join(genome_dir, "*.f*a"))
    phage_list = [".".join(op.basename(i).split(".")[:3])+"." for i in fastas]

    gff_dir = op.join(outdir, 'gff3')

    print("Writing results to gff3")
    gff_dir = write_gff3s(phage_list, gff_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_threshold=75)
    print('Done!')

if __name__ == '__main__':
    cli()
