import collections
import multiprocessing
import os
import pathlib
import shlex
import subprocess

import pandas as pd

this_directory = pathlib.Path(__file__).resolve().parent
# This should exist with source after compilation.
JUNCTION_CMD = os.path.join(this_directory, 'junction_extractor')


def get_junctions_from_bed(bed_path: str) -> dict:
    junctions = {}
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('track'):
                continue
            fields = line.split('\t')
            if len(fields) >= 6:
                chrom, start, end, name, score, strand = fields[:6]
            else:
                print(f"Error in file: {bed_path}")
                print(f"Offending line: {line}")
                raise ValueError("Invalid BED format: Expected at least 6 columns.")

            start, end = int(start), int(end)
            score = int(score)
            if start > end:
                raise Exception("Start of region cannot be greater than end of region for:\n", line)
            junctions[(chrom, start, end, strand)] = (name, score)
    return junctions

def get_junctions_multi_bed(bed_list:list, p) -> dict:
    with multiprocessing.Pool(p) as pool:
        results = pool.map(get_junctions_from_bed, bed_list)

    dd = collections.defaultdict(dict)
    for i, d in enumerate(results):
        name2 = os.path.splitext(os.path.basename(bed_list[i]))[0]
        for key, (name, score) in d.items():
            if key not in dd:
                dd[key]['samples']= set()
                dd[key]['samples'].add((name, name2, score))
                dd[key]['score'] = score
            else:
                dd[key]['samples'].add((name, name2, score))
                dd[key]['score'] = dd[key]['score'] + score

    return dd

def junction_extractor(bam_path:str, out_path:str) -> dict:
    check_for_dependency()
    name = os.path.splitext(os.path.basename(bam_path))[0]
    cmd = f"{JUNCTION_CMD} -o {out_path} {bam_path}"
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
    for _, row in df.iterrows():
        junctions[(row[0],row[1],row[2],row[5])] = (name,row[4])
    return junctions

def junction_extractor_multi_bam(bam_list:list, out_original_junctions:list, p:int) -> dict:
    with multiprocessing.Pool(p) as pool:
        results = pool.starmap(junction_extractor, zip(bam_list,out_original_junctions))

    dd = collections.defaultdict(dict)
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
    gtf = open(gtf_path, "r")

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

        if 'gene_id' not in values_dict:
            values_dict['gene_id'] = "NA"

        if 'transcript_id' not in values_dict:
            raise Exception("Exon does not contain transcript ID\n")

        transcript_id = values_dict['transcript_id']
        gene_id = values_dict['gene_id']

        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, gene_id, [[start, end]]]
        else:
            trans[transcript_id][3].append([start, end])


    for tran, [chrom, strand, gene_id, exons] in trans.items():
            exons.sort()


    junctions = collections.defaultdict(dict)
    for tran, (chrom, strand, gene_id, exons) in trans.items():
        for i in range(1, len(exons)):
            if 'transcripts' not in junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]:
                junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]['transcripts'] = [gene_id, [tran]]
            else:
                junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]['transcripts'][1].append(tran) #intron bed coordinates

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

def check_for_dependency():
    if not os.path.exists(JUNCTION_CMD):
        raise RuntimeError(f"{JUNCTION_CMD} not found.")
