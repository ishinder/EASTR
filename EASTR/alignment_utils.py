import pysam
import mappy as mp
from collections import namedtuple
import re

def get_seq(chrom,start,end,ref_fa):
    fasta = pysam.FastaFile(ref_fa)
    seq = fasta.fetch(region=chrom,start=start,end=end)
    return seq


def align_seq_pair(rseq,qseq,scoring,k,w,m):
    #TODO ambiguous bases not working
    a = mp.Aligner(seq=rseq,k=7,w=7,best_n=1,scoring=scoring,min_chain_score=m)
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

def get_alignment(chrom, jstart, jend, o5, o3, ref_fa, max_length,
                 read_length, scoring,  k, w, m):
    hits = []
    Overhang = namedtuple("Overhang","rstart rend qstart qend")
    
    #o5
    rstart = max(jstart - o5, 0)
    rend = jstart + (read_length -  o5)
    qstart = jend - o5
    qend = min(jend + (read_length -  o5), max_length)
    o5 = Overhang(rstart, rend, qstart, qend)
    
    #o3
    rstart = jend - (read_length - o3)
    rend =  min(jend + o3, max_length)
    qstart = max(jstart - (read_length - o3), 0)
    qend = jstart + o3
    o3 = Overhang(rstart, rend, qstart, qend)

    for o in (o5,o3): 
        rseq = get_seq(chrom, o.rstart, o.rend, ref_fa)
        qseq = get_seq(chrom, o.qstart, o.qend, ref_fa)
        hit = align_seq_pair(rseq, qseq, scoring,k,w,m)

        if hit:
            score = calc_alignment_score(hit, scoring,read_length)
            hits.append((hit,score))

        else:
            hits.append((None,-scoring[1]*read_length))
        
    return hits