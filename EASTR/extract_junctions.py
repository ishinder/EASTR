from collections import defaultdict, namedtuple
import pandas as pd
import os
import subprocess
import pysam
import shlex
from io import StringIO

def extract_splice_sites_bam(bam_path:str) -> dict:

    junctions = {}
    Junction = namedtuple('Junction',['score','o5','o3'])

    if not os.path.exists(bam_path + ".bai"):
        pysam.index(bam_path)
    
    cmd = "regtools junctions extract -s 0 -a 1 -m 1 -M 100000000 " + bam_path
    a = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    b = StringIO(a.communicate()[0].decode('utf-8'))

    for line in b:
        line=line.strip().split('\t')
        chrom = line[0]
        anchorStart = int(line[1])
        anchorEnd = int(line[2])
        score = int(line[4])
        strand = line[5]
        o5,o3 = [int(i) for i in line[10].split(',')] 
        juncStart = anchorStart + o5
        juncEnd = anchorEnd - o3
        junctions[(chrom,juncStart,juncEnd,strand)] = Junction(score,o5,o3)

    return junctions


def extract_splice_sites_gtf(gtf_path:str) -> dict:
    genes = defaultdict(list)
    trans = {}

    gtf =open(gtf_path, "r")

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf:
        line = line.strip()
        if line.startswith('#'):
            continue
        chrom, source, feature, start, end, score, \
                strand, frame, attributes = line.split('\t')

        start, end = int(start), int(end)

        if start > end:
            raise Exception("Start of region can not be greater than end of region for:\n",line)

        if feature != 'exon':
            continue

        values_dict = {}
        for attr in attributes.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            raise Exception("Exon does not contain transcript or gene ID\n")

        transcript_id = values_dict['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[start, end]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            trans[transcript_id][2].append([start, end])


    for tran, [chrom, strand, exons] in trans.items():
            exons.sort()


    junctions = defaultdict(list)
    for tran, (chrom, strand, exons) in trans.items():
        for i in range(1, len(exons)):
            junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)].append(tran) #intron bed coordinates
    
    return junctions

if __name__ == '__main__':
    bam_path='tests/data/ERR188044_chrX.bam'
    gtf_path = '/ccb/salz8-2/shinder/projects/EASTR/EASTR/tests/data/chrX.gtf'
    junctions_bam = extract_splice_sites_bam(bam_path)
    junctions_gtf = extract_splice_sites_gtf(gtf_path)
