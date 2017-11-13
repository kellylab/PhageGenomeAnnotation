from __future__ import print_function
import click
import os.path as op
import os
import glob

from phage_ublast import run_ublasts
from phage_prodigal import run_prodigals
from phage_trnascan import run_trna_scans, prep_outdir
from phage_crt import run_crts
from nvp_add_ips import merge_gff3_and_ips
from nvp_output_scripts import write_gff3_file

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    '''functions to output gff3s'''
    pass

def write_gff3s(phage_list, output_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, prefix="NVP", cov_threshold=75, overwrite=False):
    prep_outdir(output_dir)
    for phage in phage_list:
        output_file = op.join(output_dir, '{}.gff3'.format(phage))
        if op.exists(output_file) and overwrite == False:
            return output_dir

        # write_gff3_file(phage, output_file, prod_path, faa_path, genome_path, blast_path, trna_path, crt_path, cov_thresh=75):
        else:
            write_gff3_file(phage, output_file, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_thresh=75, prefix=prefix)
    return output_dir


@cli.command("phage-list-runall", short_help="provide space separated list of phages to blast")
@click.argument('phage-list', nargs=-1)
@click.option('--outdir', help="where to send blast outputs")
@click.option('--genome-dir', help="where to find genomic contigs in fasta format")
@click.option('--blast_databasedir', help="where to find the blast databases", default='/nobackup1/jbrown/annotation/databases', show_default=True)
@click.option('--ublast_path', help="where ublast executable found", default='/home/sbiller/usearch7.0.1090_i86linux64', show_default=True)
@click.option('--ublast_evalue', help='evalue to use for ublast', default='1e-5', show_default=True)
@click.option('--path-to-crt', help='location of crt executable', default="/home/jbrown/programs/CRT1.2-CLI.jar", show_default=True)
@click.option("--nvp",
                help="True if this phage is formatted like a nahant vibriophage (#.###.X) e.g. 1.028.A, otherwise, phage name assumed to be everything before first '.' e.g. for file /genomedir/phage1.fasta, the name is 'phage1'",
                default=True, show_default=True
                )
def run_all(phage_list, genome_dir, outdir, blast_databasedir, ublast_path, ublast_evalue, path_to_crt):
    print("Running prodigal now")

    if nvp is True:
        prefix = "NVP"
    else:
        prefix = "CDS"
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
    gff_dir = write_gff3s(phage_list, gff_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_threshold=75, prefix=prefix)


@cli.command("genome-dir-runall", short_help="indicate directory where all genomes are found as fasta files (all files in directory will be annotated)")
@click.option('--outdir', help="where to send outputs")
@click.option('--genome-dir', help="where to find genomic contigs in fasta format")
@click.option('--blast_databasedir', help="where to find the blast databases", default='/nobackup1/jbrown/annotation/databases', show_default=True)
@click.option('--ublast_path', help="where ublast executable found", default='/home/sbiller/usearch7.0.1090_i86linux64', show_default=True)
@click.option('--ublast_evalue', help='evalue to use for ublast', default='1e-5', show_default=True)
@click.option('--path-to-crt', help='location of crt executable', default="/home/jbrown/programs/CRT1.2-CLI.jar", show_default=True)
@click.option("--nvp",
                help="True if this phage is formatted like a nahant vibriophage (#.###.X) e.g. 1.028.A, otherwise, phage name assumed to be everything before first '.' e.g. for file /genomedir/phage1.fasta, the name is 'phage1'",
                default=True, show_default=True
                )
def run_all_dir(genome_dir, outdir, blast_databasedir, ublast_path, ublast_evalue, path_to_crt, nvp):
    fastas = glob.glob(op.join(genome_dir, "*.f*a"))

    if nvp is True:
        phage_list = [".".join(op.basename(i).split(".")[:3]) for i in fastas]
        prefix = "NVP"
    else:
        phage_list = [i.split(".")[0] for i in fastas]
        prefix = "CDS"

    print("Running prodigal now")
    #run_prodigals(phage_names, outdir, genome_dir)
    prod_dir, prot_dir, fna_dir = run_prodigals(phage_list, outdir, genome_dir)

    blast_dir = op.join(outdir, "blasts")
    print("Running blasts now")
    #run_ublasts(phage_list, outdir, databasedir, proteindir, ublast_path='/home/sbiller/usearch7.0.1090_i86linux64', evalue='1e-5'):
    blast_dir = run_ublasts(phage_list, blast_dir, blast_databasedir, prot_dir, ublast_path, ublast_evalue)

    trna_dir = op.join(outdir, 'trna')
    print("Running trna scan now")
    trna_dir = run_trna_scans(phage_list, trna_dir, genome_dir)

    crt_dir = op.join(outdir, 'crt')
    print("Looking for CRISPRs using crt")
    crt_dir = run_crts(phage_list, crt_dir, genome_dir, path_to_crt=path_to_crt)

    gff_dir = op.join(outdir, 'gff3')
    print("Writing results to gff3")
    gff_dir = write_gff3s(phage_list, gff_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_threshold=75, prefix=prefix)


@cli.command("write-gff3s", short_help='write gff3 files given all appropriate input directories')
@click.option("--genome-dir", help="location of genomic contigs in fasta format")
@click.option("--outdir", help="parent directory of all other results directories, assumes that trnascan, crt and blasts have already been run and have default naming structure")
@click.option("--ips-dir", help="interproscan results location", default=None)
@click.option("--blast-dir", help="location of blast outputs for annotation", default=None)
@click.option("--trna-dir", help="location of trnascan results if not in outdir", default=None)
@click.option("--crt-dir", help="location of crt results if not in outdir", default=None)
@click.option("--prod-dir", help="location of prodigal .gff files if not in outdir", default=None)
@click.option("--prot-dir", help="location of translated orfs if not in outdir", default=None)
@click.option("--overwrite", help="Ignore existing gff3 files if True", default=False)
@click.option("--nvp",
                help="True if this phage is formatted like a nahant vibriophage (#.###.X) e.g. 1.028.A, otherwise, phage name assumed to be everything before first '.' e.g. for file /genomedir/phage1.fasta, the name is 'phage1'",
                default=True)
def write_from_genomedir(genome_dir, outdir, ips_dir, blast_dir, trna_dir, crt_dir, prod_dir, prot_dir, overwrite, nvp):
    '''
    nvp_annotations.py write-gff3s --genome-dir /nobackup1/jbrown/nvp_for_ncbi/genomes/ --outdir /nobackup1/jbrown/nvp_for_ncbi/ \
    --ips-dir /nobackup1/jbrown/annotation/protein_ips_calls/ --blast-dir /nobackup1/jbrown/annotation/blasts/ \
    --trna-dir /nobackup1/jbrown/annotation/trna --crt-dir /nobackup1/jbrown/annotation/crt/ \
    --prod-dir /nobackup1/jbrown/annotation/genes/ --prot-dir /nobackup1/jbrown/annotation/proteins/
    '''

    if blast_dir is None: blast_dir = op.join(outdir, "blasts")
    if trna_dir is None: trna_dir = op.join(outdir, 'trna')
    if crt_dir is None: crt_dir = op.join(outdir, 'crt')
    if prod_dir is None: prod_dir = op.join(outdir, 'genes')
    if prot_dir is None: prot_dir = op.join(outdir, 'proteins')

    for i in [blast_dir, trna_dir, crt_dir, prod_dir, prot_dir]: assert op.exists(i), "Please make sure that {i} exists with a result per genome".format(i=i)

    fastas = glob.glob(op.join(genome_dir, "*.f*a"))
    if nvp is True:
        phage_list = [".".join(op.basename(i).split(".")[:3]) for i in fastas]
        prefix = "NVP"
    else:
        phage_list = [i.split(".")[0] for i in fastas]
        prefix = "CDS"

    gff_dir = op.join(outdir, 'gff3')

    print("Writing results to gff3")
    gff_dir = write_gff3s(phage_list, gff_dir, prod_dir, prot_dir, genome_dir, blast_dir, trna_dir, crt_dir, cov_threshold=75, overwrite=overwrite, prefix=prefix)

    if ips_dir is not None:
        ips_merge_out = op.join(outdir, 'gff_ips_added')
        prep_outdir(ips_merge_out)
        merge_gff3_and_ips(ips_dir, gff_dir, ips_merge_out)

if __name__ == '__main__':
    cli()
