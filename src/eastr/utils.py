import collections
import os

import pysam


def index_fasta(ref_fa):
    if not os.path.exists(f"{ref_fa}.fai"):
        pysam.faidx(ref_fa)

#Make a new directory
def make_dir(path):
    directory = os.path.join(path)
    os.makedirs(directory,exist_ok=True)


def get_chroms_list_from_fasta(ref_fa):
    fasta=pysam.FastaFile(ref_fa)
    chroms = list(fasta.references)
    chrom_sizes = collections.defaultdict(int)
    for chrom in chroms:
        chrom_sizes[chrom] = fasta.get_reference_length(chrom)
    return chrom_sizes

def check_directory_or_file(path:str) -> str:
    if os.path.splitext(os.path.basename(path))[1]!='':
        return 'file'
    else:
        return 'dir'
