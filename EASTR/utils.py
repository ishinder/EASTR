import collections
import os
import pathlib

import pkg_resources
import pysam


def check_cram_bam(bam_path: str):
    path_parameter = bam_path.split(".")[-1]
    if path_parameter in {"bam", "cram"}:
        return path_parameter
    else:
        raise Exception("Input must be a cram or a bam file")
    

def index_fasta(ref_fa):
    if not os.path.exists(f"{ref_fa}.fai"):
        pysam.faidx(ref_fa)

# UNUSED
def index_bam(bam_path):
    if check_cram_bam(bam_path) == "bam":
        if not os.path.exists(bam_path + ".bai"):
            pysam.index(bam_path)
    else:
        if not os.path.exists(bam_path + ".crai"):
            pysam.index(bam_path)


# UNUSED
def get_read_length_from_bam(bam_file):
    read_length = read_sum = 0
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for bam_length in samfile.head(1000000):
        read_length += 1
        read_sum += bam_length.query_length
    if read_length == 0:
        return 0
    return int(read_sum / read_length)

def make_dir(path):
    d = pathlib.Path(path)
    d.mkdir(exist_ok=True)

# UNUSED
def get_chroms_list_from_bam(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    chroms = list(samfile.references)
    chrom_sizes = collections.defaultdict(int)
    for chrom in chroms:
        chrom_sizes[chrom] = samfile.get_reference_length(chrom)
    return chrom_sizes


def get_chroms_list_from_fasta(ref_fa):
    fasta = pysam.FastaFile(ref_fa)
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

def check_directory_or_file(path: str) -> str:
    p = pathlib.Path(path)
    if p.is_dir():
        return "file"
    elif p.is_file():
        return "dir"
    else:
        return ""

def is_directory(path: str) -> bool:
    p = pathlib.Path(path)
    return p.is_dir()

def is_file(path: str) -> bool:
    p = pathlib.Path(path)
    return p.is_file()

def get_bam_list(bam_file):
    # If a single bam file is provided
    extension = os.path.splitext(os.path.basename(bam_file))[1]
    if extension in {".bam",".cram",".sam"}:
        bam_list = [bam_file]
    else:
        with open(bam_file) as opened_bam_file:
            bam_list = []
            for line in opened_bam_file:
                bam_file = line.rstrip()
                if os.path.isfile(bam_file):
                    bam_list.append(bam_file)
                else:
                    raise ValueError("Input must be a bam file or a file "
                                         "containing a list of bam files")                   

    return bam_list

def get_bed_list(bed_file):
    extension = os.path.splitext(os.path.basename(bed_file))[1]
    if extension in {".bed"}:
        bed_list = [bed_file]
    else:
        with open(bed_file) as opened_bed_file:
            bed_list = []
            for line in opened_bed_file:
                bed_file = line.rstrip()
                if os.path.isfile(bed_file):
                    bed_list.append(bed_file)
                else:
                    raise ValueError("Input must be a bed file or a file "
                        "containing a list of bed files")
    return bed_list

def get_minimap_scoring_list(args) -> list:
    gap_open_penalty = args.O
    gap_ext_penalty = args.E
    mismatch_penalty = args.B
    match_score = args.A
    ambiguous_score = args.scoreN

    scoring = [
        match_score,
        mismatch_penalty,
        gap_open_penalty[0],
        gap_ext_penalty[0],
        gap_open_penalty[1],
        gap_ext_penalty[1],
        ambiguous_score,
        ]

    return scoring
