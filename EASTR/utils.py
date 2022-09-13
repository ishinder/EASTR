from subprocess import run
import os
import shlex
import pysam

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

