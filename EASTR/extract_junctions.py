from collections import defaultdict, namedtuple
import collections
import pickle
from posixpath import basename
import tempfile
import pandas as pd
import os
import subprocess
import pysam
import shlex
from io import StringIO


class Alignment:
    def __init__(self,alignment):
        self.query_name = alignment.query_name
        self.reference_name = alignment.reference_name
        self.reference_start = alignment.reference_start
        self.NH = alignment.get_tag("NH")
        self.is_proper_pair = alignment.is_proper_pair
        self.is_read1 = alignment.is_read1
        self.next_reference_start = alignment.next_reference_start
        self.cigarstring = alignment.cigarstring
        self.cigar = alignment.cigar


def extract_splice_sites_bam_regtools(bam_path:str) -> dict:
    junctions = {}
    # Junction = namedtuple('Junction',['score','o5','o3'])

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
        junctions[(chrom,juncStart,juncEnd,strand)] = (score,o5,o3)

    return junctions

def extract_splice_sites_bam_custom(bam_path,tmp_path=os.getcwd()):

    prefix = basename(bam_path).split('.')[0]
    suffix = ".junctions.pickle"
    

    samfile = pysam.AlignmentFile(bam_path)
    junctions=collections.OrderedDict() 

    for alignment in samfile.fetch(until_eof=True, multiple_iterators=True):
        currentloc = alignment.pos
        for i,cigarop in enumerate(alignment.cigar):
            if (cigarop[0]==4): #substitution or insertion in query
                continue
            if (cigarop[0]==1):
                continue
            if(cigarop[0]==3):
                key = (alignment.reference_name,currentloc,currentloc+cigarop[1],alignment.get_tag('XS'))
                a = Alignment(alignment)
                o5 = a.cigar[i-1][1]
                o3 = a.cigar[i+1][1]
                if key not in junctions:
                    junctions[key]=[(a, o5, o3)]
                else:
                    junctions[key].append((a, o5, o3))         
            currentloc=currentloc+cigarop[1]
    
    tmpfile = tempfile.NamedTemporaryFile(dir=tmp_path, prefix=prefix, suffix=suffix, delete=False)
    pickle.dump(junctions, tmpfile)
    tmpfile.close()
    
    return junctions



def extract_splice_sites_gtf(gtf_path:str) -> dict:
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
