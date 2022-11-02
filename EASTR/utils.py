from subprocess import run
import os
import shlex
from xml.sax.handler import feature_namespaces
import pysam
from posixpath import dirname
import collections

def index_fasta(ref_fa):
    if not os.path.exists(f"{ref_fa}.fai"):
        pysam.faidx(ref_fa)

def index_bam(bam):
    if not os.path.exists(f"{bam}.bai"):
        pysam.index(bam)

def get_read_length_from_bam(bam):
    read_lengths = []
    samfile = pysam.AlignmentFile(bam, "rb")
    for l in samfile.head(1000000):
        read_lengths.append(l.query_length)
        
    return int(sum(read_lengths) / len(read_lengths))

#Make a new directory
def make_dir(path):
    directory = dirname(path)
    isExist = os.path.exists(directory)   
    if not isExist:
        os.makedirs(directory)

def get_chroms_list_from_bam(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    chroms = list(samfile.references)
    chrom_sizes = collections.defaultdict(int)
    for chrom in chroms:
        chrom_sizes[chrom] = samfile.get_reference_length(chrom)
    return chrom_sizes


def get_chroms_list_from_fasta(ref_fa):
    fasta=pysam.FastaFile(ref_fa)
    chroms = list(fasta.references)
    chrom_sizes = collections.defaultdict(int)
    for chrom in chroms:
        chrom_sizes[chrom] = fasta.get_reference_length(chrom)
    return chrom_sizes
