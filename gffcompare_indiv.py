import glob
import shlex
import subprocess

from collections import defaultdict
from EASTR import utils
import glob
import argparse
import os

def parse_args(arglist):
    parser = argparse.ArgumentParser(
    prog="EASTR tests: individual GFF compare",
    description="Run gffcompare for original and filtered BAM samples")

    parser.add_argument("--dir", help="base directory", default=os.getcwd())
    parser.add_argument("-p", help="number of processes", type=int)
    parser.add_argument("-r", help="reference GTF for sensitivity/precision calculation")
    parser.add_argument("--outdir", help="base directory", default=f"{os.getcwd()}/gffcompare/")
    parser.add_argument("--filters", help="keywords to filter by, comma seperated list", default=None)
    parser.add_argument("--glob_by", help="regular expression string to select specific files", default=None)

    return parser.parse_args()

def print_gffcompare_results(filter_by, outdir, glob_by):

    cmd=f"grep \"{filter_by}\" {outdir}/*/{glob_by}stats | awk -F'original/|filtered/' '{{print $2}}' | sort -nk1,1 | column -t"
    os.system(cmd)
    # grep "Transcript level" /ccb/salz2/shinder/projects/EASTR_tests/maize/gffcompare//*/*stats | awk -F'original/|filtered/' '{print $2}' | sort -nk1,1 | column -t


def main(arglist=None):
    args = parse_args(arglist)
    p = args.p
    BASEDIR = args.dir
    ref_gtf = args.r
    outdir = args.outdir
    glob_by = args.glob_by
    filters = args.filters

    if glob_by is None:
        glob_by= "*"

    if filters is None:
        filters="Transcript level, Intron level, Missed exons, Novel introns, Novel exons, Novel loci, Matching transcripts"

    filters = filters.split(',')

    BASEDIR="/ccb/salz8-3/shinder/projects/EASTR_tests/CD4"
    outdir=f"{BASEDIR}/gffcompare"
    ref_gtf="/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf"

    gtf_list = {}
    for t in ["original","filtered"]:
        odir = f"{outdir}/{t}"
        utils.make_dir(odir)
        os.chdir(odir)
        gtf_list[t] = glob.glob(f"{BASEDIR}/GTF/{t}/*{glob_by}*.gtf")
        for file in gtf_list[t]:

            name=os.path.basename(file).split('.')[0]
            if not os.path.exists(f"{odir}/{name}.stats"):
                cmd = f"gffcompare -r {ref_gtf} {file} -o {name}"
                subprocess.run(shlex.split(cmd))

    # filters=["Novel loci"]
    for filter_by in filters:
        print_gffcompare_results(filter_by, outdir, glob_by)


if __name__ == '__main__':
    main()

