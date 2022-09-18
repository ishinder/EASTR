
from typing import Dict
import pandas as pd
import pysam
from multiprocessing import Pool, Manager
import mappy as mp
import collections
import os
import subprocess
import shlex
from io import StringIO
import re
from EASTR import utils
import numpy as np
from multiprocessing import Pool
from collections import namedtuple

# def get_introns_from_bam(samfile):
#     introns={}
#     for read in samfile.fetch(until_eof=True):
#         currentloc = read.pos
#         for i,cigarop in enumerate(read.cigar):
#             if (cigarop[0]==4): #substitution or insertion in query
#                 continue
#             if (cigarop[0]==1):
#                 continue
#             if(cigarop[0]==3):
#                 key = (samfile.getrname(read.tid),currentloc,currentloc+cigarop[1])
#                 o5 = read.cigar[i-1][1]           
#                 o3 = read.cigar[i+1][1]  
#                 if key not in introns:
#                     introns[key]=[(read,o5,o3)]
#                 else:
#                     introns[key].append((read,o5,o3))            
#             currentloc=currentloc+cigarop[1]
#     return introns

#TODO: regtools is 2.2X faster than the function above (80 vs. 173 seconds).
def get_introns_regtools(bamfile:str, chrom=None) -> pd.DataFrame:
    if not os.path.exists(bamfile + ".bai"):
        pysam.index(bamfile)
    
    if chrom is None:
        cmd = "regtools junctions extract -s 0 -a 1 -m 1 -M 100000000 " + bamfile
        
    else: 
        cmd = f"regtools junctions extract -s 0 -a 1 -m 1 -M 100000000 -r {chrom}" + bamfile
        
    a = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    b = StringIO(a.communicate()[0].decode('utf-8'))

    colnames = ['chrom','anchorStart','anchorEnd','name','score','strand','anchorSize']
    df = pd.read_csv(b, header=None, index_col=None, sep='\t', usecols=[0,1,2,3,4,5,10],names=colnames)
    
    df['5o'] = df['anchorSize'].str.split(',').str[0].astype(int)
    df['3o'] = df['anchorSize'].str.split(',').str[1].astype(int)
    df['juncStart'] = df['anchorStart'] + df['5o']
    df['juncEnd'] = df['anchorEnd'] - df['3o']
    
    df = df.set_index(["chrom","juncStart","juncEnd"]).sort_index()
    df.drop(labels = 'anchorSize',axis=1,inplace=True)

    return df


def get_seq(chrom,start,end,ref_fa):
    fasta = pysam.FastaFile(ref_fa)
    seq = fasta.fetch(region=chrom,start=start,end=end)
    return seq

def align_seq_pair(rseq,qseq,scoring,k,w):
    #TODO ambiguous bases not working
    a = mp.Aligner(seq=rseq,k=7,w=7,best_n=1,scoring=scoring)
    itr = list(a.map(qseq,MD=True,cs=True))
    
    if not itr:
        return
    
    for hit in itr:

        if hit.strand != 1:
            continue

        if not hit.is_primary:
            continue

        else:
            return hit #returns the first primary hit


def calc_alignment_score(hit,scoring,read_length):
    #TODO verify AS/alignment store calc
    matches = hit.mlen
    gap_penalty = 0
    cs = hit.cs
    
    #gaps
    p = re.compile('[\\-\\+]([atgc]+)')
    m = p.findall(cs)
    gaps = len(m)
    for gap in m:
        gap_len = len(gap)
        gap_penalty += min(scoring[2] + (gap_len - 1) * scoring[3],
                            scoring[4] + (gap_len - 1) * scoring[5])  
        
    #mismatches
    p = re.compile('\\*([atgc]+)')
    m = p.findall(cs)
    mismatches = len(m)
        
    AS = matches*scoring[0] - (mismatches)*scoring[1] - gap_penalty

    return AS


# def get_alignment(chrom, start, end, rstart, rend, qstart, qend, ref_fa, scoring):
#     rseq = get_seq(chrom, rstart, rend, ref_fa)
#     qseq = get_seq(chrom, qstart, qend, ref_fa)
#     hit = align_seq(rseq,qseq)

#     if hit:                
#         row = [ chrom,
#                 start,
#                 end,
#                 rstart + hit.r_st, 
#                 rstart + hit.r_en,
#                 qstart + hit.q_st,
#                 qstart + hit.q_en,
#                 hit.mlen,
#                 calc_alignment_score(hit, scoring)]

#         return row
#     else:
#         return None

def get_alignment(chrom, jstart, jend, o5, o3, scoring, ref_fa, read_length,k,w):
    hits = []
    Overhang = namedtuple("Overhang","rstart rend qstart qend")
    
    #o5
    rstart = jstart - o5
    rend = jstart + (read_length -  o5)
    qstart = jend - o5
    qend = jend + (read_length -  o5)
    o5 = Overhang(rstart, rend, qstart, qend)
    
    #o3
    rstart = jend - (read_length - o3)
    rend =  jend + o3
    qstart = jstart - (read_length - o3)
    qend = jstart + o3
    o3 = Overhang(rstart, rend, qstart, qend)
    

    for o in (o5,o3): 
        rseq = get_seq(chrom, o.rstart, o.rend, ref_fa)
        qseq = get_seq(chrom, o.qstart, o.qend, ref_fa)
        hit = align_seq_pair(rseq, qseq, scoring,k,w)

        if hit:
            score = calc_alignment_score(hit, scoring,read_length)
            hits.append((hit,score))

        else:
            hits.append((None,-4*read_length))
        
    return hits

    

def run_junctions(bam, scoring, ref_fa, read_length,k,w):
    introns = get_introns_regtools(bam)
    
    for junc,row in introns.iterrows():
        chrom = junc[0]
        jstart = junc[1]
        jend = junc[2]
        o5 = row["5o"]
        o3 = row["3o"]
        hits = get_alignment(chrom, jstart, jend, o5, o3, scoring, ref_fa, read_length,k,w)
        
    
        introns.loc[junc,"5o_score"] = hits[0][1]
        introns.loc[junc,"3o_score"] = hits[1][1]
        introns.loc[junc,"max_o_score"] = np.nanmax([hits[0][1],hits[1][1]])
        
    return introns

#def run_junctions_parallel(bam,scoring,ref_fa,read_length,k,w)
        
    

if __name__ == '__main__':
    import time
    start = time.time()
    
    bam = "tests/data/ERR188044_chrX.bam"
    ref_fa = "tests/data/chrX.fa"
    scoring=[2,4,4,2,24,1,1]
    read_length = utils.get_read_length_from_bam(bam)
    k = 7
    w = 7
    #chroms = [x['SN'] for x in samfile.header['SQ']]

    introns = run_junctions(bam, scoring, ref_fa, read_length,k,w)
    
    end = time.time()
    print(f"took {end-start} seconds")
    
    introns[~introns['score'].isna()]