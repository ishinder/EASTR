import collections
import os

import pysam


def index_fasta(ref_fa):
    if not os.path.exists(f"{ref_fa}.fai"):
        pysam.faidx(ref_fa)

def make_dir(path: str):
    """Make a directory, if it exists just skip."""
    directory = os.path.join(path)
    os.makedirs(directory, exist_ok=True)


def get_chroms_list_from_fasta(ref_fa):
    fasta=pysam.FastaFile(ref_fa)
    chroms = list(fasta.references)
    chrom_sizes = collections.defaultdict(int)
    for chrom in chroms:
        chrom_sizes[chrom] = fasta.get_reference_length(chrom)
    return chrom_sizes

def sanitize_name(path: str) -> str:
    """Return base name from path without extension.
    Example:
        input: /var/tmp/file.ext
        output: file
    """
    base_name = os.path.basename(path)
    name_without_extension, _ = os.path.splitext(base_name)
    return name_without_extension

def get_file_extension(path: str) -> str:
    """Return file extension name from a file path.
    Example:
        input: /var/tmp/file.ext
        output: ext
    """
    base_name = os.path.basename(path)
    _, extension = os.path.splitext(base_name)
    return extension

def sanitize_and_update_extension(path, extension) -> str:
    base_name = sanitize_name(path)
    return f"{base_name}{extension}"
