from __future__ import print_function
import glob
import os
import os.path as op
import pandas as pd
import click

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    '''merge gff and ips'''
    pass

def ips_dict(ips):
    ips_columns = ["ProteinAcc",
               "SeqMD5",
               "SeqLen",
               "Analysis",
               "SignatureAccession",
               "SignatureDescription",
               "StartLoc",
               "StopLoc",
               "Score",
               "Status",
               "Date",
               "InterProAnnotationAcc",
               "InterproAnnotationDesc",
               "GOAnnotation"]
    df = pd.read_csv(ips, sep="\t", index_col=False, names=ips_columns)
    best_hits = df[df.groupby(['ProteinAcc'])['Score'].transform(min)==df['Score']].drop_duplicates()
    ips_anns = {}

    for i, l in best_hits.iterrows():
        pid = l.ProteinAcc
        source = "InterPro"
        if pd.isnull(l.InterproAnnotationDesc):
            if l.Analysis == "TMHMM" or l.Analysis == "Coils":
                desc = "{element} containing protein".format(element=l.SignatureAccession)
            else:
                if pd.isnull(l.SignatureDescription):
                    continue
                else:
                    desc = l.SignatureDescription
            ips_anns[pid] = [desc, source]
        else:
            desc = l.InterproAnnotationDesc
            ipsid = l.InterProAnnotationAcc
            ips_anns[pid] = [desc,source, ipsid]

    return ips_anns

class GffLine():
    def __init__(self, line):
        vec = [i.rstrip() for i in line.strip().split("\t")]
        self.name = vec[0]
        self.method = vec[1]
        self.etype = vec[2]
        self.cstart = vec[3]
        self.cstop = vec[4]
        self.dot = vec[5]
        self.strand = vec[6]
        self.something = vec[7]
        self.notes = ";".join([i for i in vec[8].split(";") if "ID=" not in i and "Name=" not in i])
        self.pid = [i for i in vec[8].split(";") if "ID=" in i][0]
        self.desc = [i for i in vec[8].split(";") if "Name=" in i][0]
        self.key = "{name}_{number}".format(name=self.name, number=int(self.pid.split("_")[-1]))

    def construct_note(self):
        notes = ";".join([self.pid, self.desc, self.notes])

    def print_line(self):
        notes = ";".join([self.pid, self.desc, self.notes])
        line = "\t".join([self.name, self.method, self.etype, self.cstart,
                          self.cstop, self.dot, self.strand, self.something, notes])
        return line

    def change_id(self, newid):
        self.pid = "ID={newid}".format(newid=newid)

    def change_desc(self, newdesc):
        self.desc = "Name={newdesc}".format(newdesc=newdesc)


def find_file_matches(ipslist, gfflist):
    for i in ipslist:
        phage = ".".join(op.basename(i).split(".")[:3])
        try:
            ann = [j for j in gfflist if phage in j][0]
        except:
            print("gff for {} not found".format(phage))
            ann = None
        yield phage, i, ann


def combine_ips_gff3(ips, gff3, outfile):
    ipsdict = ips_dict(ips)
    with open(gff3) as infile, open(outfile, "w") as oh:
        for l in infile:
            if "tRNAScan" not in l and "crispr" not in l.lower():
                line = GffLine(l)
                if line.key in ipsdict.keys():
                    if "hypothetical" in line.desc or "unknown" in line.desc:
                        line.change_desc(ipsdict[line.key][0])
                    if len(ipsdict[line.key]) == 3:    
                        line.notes += '; note="InterPro:{ips}"'.format(ips=ipsdict[line.key][-1])
                print(line.print_line(), file=oh)
            else:
                print(l, file=oh)

def merge_gff3_and_ips(ips_dir, gff3_dir, outdir):
    if op.exists(outdir) == False:
        os.mkdir(outdir)

    ips_list = glob.glob(op.join(ips_dir, "*.tsv"))
    gff_list = glob.glob(op.join(gff3_dir, "*"))

    for phage, ips, gff in find_file_matches(ips_list, gff_list):
        if gff is None:
            print("WARNING: could not find gff file for {}, skipping pairing".format(phage))
            continue
        out_gff = op.join(outdir, "{}_ips.gff3".format(phage))
        print("merging {} and {} into {} for phage {}".format(ips, gff, out_gff, phage))
        combine_ips_gff3(ips, gff, out_gff)
    print("done!")

@cli.command('merge-gff-ips-files', short_help='merge a gff3 file and a ips tsv')
@click.argument('gff')
@click.argument('ips')
@click.argument('outfile')
def merge_gff_ips_files(gff, ips, outfile):
    combine_ips_gff3(ips, gff, outfile)


@cli.command("merge-dirs", short_help="merge gff and ips tsv outputs for nvp collection")
@click.option('--ips-dir', help="directory containing interproscan tsv results")
@click.option('--gff3-dir', help="directory containing gff3 annotations")
@click.option('--outdir', help="directory containing combined gff3 outputs", default="./gff3_ips")
def merge_dirs(ips_dir, gff3_dir, outdir):
    merge_gff3_and_ips(ips_dir, gff3_dir, outdir)

if __name__ == '__main__':
    cli()
