from subprocess import run
import os
import shlex
from xml.sax.handler import feature_namespaces

import pkg_resources
import pysam
from posixpath import dirname
import collections

def check_cram_bam(bam_path):
    if bam_path.split('.')[-1]=="bam":
        return 'bam'
    if bam_path.split('.')[-1]=="cram":
        return 'cram'
    else:
        raise Exception("Input must be a cram or a bam file")
    

def index_fasta(ref_fa):
    if not os.path.exists(f"{ref_fa}.fai"):
        pysam.faidx(ref_fa)

def index_bam(bam_path):
    if check_cram_bam(bam_path) == 'bam':
        if not os.path.exists(bam_path + ".bai"):
            pysam.index(bam_path)

    else :
        if not os.path.exists(bam_path + ".crai"):
            pysam.index(bam_path)


def get_read_length_from_bam(bam):
    read_lengths = []
    samfile = pysam.AlignmentFile(bam, "rb")
    for l in samfile.head(1000000):
        read_lengths.append(l.query_length)
        
    return int(sum(read_lengths) / len(read_lengths))

#Make a new directory
def make_dir(path):
    directory = os.path.join(path)
    os.makedirs(directory,exist_ok=True)

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


def get_package_path():
    eastr_package_name = "EASTR"
    eastr_path = pkg_resources.resource_filename(eastr_package_name, "")
    parent_dir = os.path.dirname(eastr_path)
    return parent_dir

def get_vacuum_path():
    eastr_path = get_package_path()
    vacuum_path = os.path.join(eastr_path, "utils", "vacuum")
    return vacuum_path

def get_junction_extractor_path():
    eastr_path = get_package_path()
    junction_extractor_path = os.path.join(eastr_path, "utils", "junction_extractor")
    return junction_extractor_path

def check_directory_or_file(path:str) -> str:
    if os.path.splitext(os.path.basename(path))[1]!='':
        return 'file'
    else:
        return 'dir'

