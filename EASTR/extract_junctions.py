from collections import defaultdict, namedtuple
import collections
import multiprocessing
import pickle
from posixpath import basename
import tempfile
import pandas as pd
import os
import subprocess
import pysam
import shlex
from io import StringIO
from EASTR import utils



def get_junctions_from_bed(bed_path:str) -> dict:
    junctions = {}
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('track'):
                continue
            chrom, start, end, name, score, strand = line.split('\t')
            start, end = int(start), int(end)
            if start > end:
                raise Exception("Start of region can not be greater than end of region for:\n",line)
            junctions[(chrom, start, end, strand)] = (name, score)
    return junctions

def junction_extractor(bam_path:str, out_path:str) -> dict:
    name = os.path.splitext(os.path.basename(bam_path))[0]
    je = utils.get_junction_extractor_path()
    cmd = je + " -o " + out_path + " " + bam_path
    a = subprocess.Popen(shlex.split(cmd), stdout=subprocess.DEVNULL)
    b = a.communicate() 

    with open(out_path, 'r') as f:
        first_line = f.readline().strip()
    if 'track name=junctions' in first_line:
        skip = 1
    else:
        skip = 0

    df = pd.read_csv(out_path, sep='\t', header=None, skiprows=skip, comment='#')
    
    junctions = {}
    for index, row in df.iterrows():
        junctions[(row[0],row[1],row[2],row[5])] = (name,row[4])
    return junctions

def junction_extractor_multi_bam(bam_list:list, out_original_junctions:list, p:int) -> dict:
    pool = multiprocessing.Pool(p)
    results = pool.starmap(junction_extractor, zip(bam_list,out_original_junctions))
    pool.close
    dd = defaultdict(dict)

    for d in results:
        for key, value in d.items():
            if key not in dd:
                dd[key]['samples']= set()
                dd[key]['samples'].add(value)
                dd[key]['score'] = value[1]
            else:
                dd[key]['samples'].add(value)
                dd[key]['score'] = dd[key]['score'] + value[1]

    return dd

def extract_splice_sites_gtf(gtf_path:str) -> dict:
    trans = {}
    gtf =open(gtf_path, "r")

    for line in gtf:
        line = line.strip()
        if line.startswith('#'):
            continue
        chrom, source, feature, start, end, score, \
                strand, frame, attributes = line.split('\t')

        start, end = int(start), int(end)

        if feature != 'exon':
            continue
        
        if start > end:
            raise Exception("Start of region can not be greater than end of region for:\n",line)

        values_dict = {}
        for attr in attributes.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'transcript_id' not in values_dict:
            raise Exception("Exon does not contain transcript ID\n")

        transcript_id = values_dict['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[start, end]]]
        else:
            trans[transcript_id][2].append([start, end])


    for tran, [chrom, strand, exons] in trans.items():
            exons.sort()


    junctions = defaultdict(dict)
    for tran, (chrom, strand, exons) in trans.items():
        for i in range(1, len(exons)):
            if 'transcripts' not in junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]:
                junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]['transcripts'] = [tran]
            else:
                junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]['transcripts'].append(tran) #intron bed coordinates
    
    return junctions

# if __name__ == '__main__':
#     bam_path='tests/data/ERR188044_chrX.bam'
#     gtf_path = '/ccb/salz8-2/shinder/projects/EASTR/EASTR/tests/data/chrX.gtf'
#     junctions_bam = extract_splice_sites_bam_regtools(bam_path)
#     junctions_gtf = extract_splice_sites_gtf(gtf_path)

#     import time
#     bam_path='/ccb/salz8-2/shinder/projects/EASTR_tests/chess_brain/BAM/original/R2947_C48YFACXX.bam'
#     start = time.time()
#     j = extract_splice_sites_bam_custom(bam_path)
#     end = time.time()
#     print(f"took {(end-start)/60} mins to get junctions")

#     start = time.time()
#     x=extract_splice_sites_bam_regtools(bam_path)
#     end = time.time()
#     print(f"took {(end-start)/60} mins to get junctions")
