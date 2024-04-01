import os
import re
import shlex
import subprocess
import tempfile

import mappy as mp

def get_seq(chrom,start,end,pysam_fa):
    seq = pysam_fa.fetch(region=chrom,start=start,end=end)
    return seq


def align_seq_pair(rseq:str, qseq:str, scoring:list, k:int, w:int, m:int, best_n=1):
    #TODO ambiguous bases not working
    a = mp.Aligner(seq=rseq,k=k,w=w,best_n=best_n,scoring=scoring,min_chain_score=m)
    itr = list(a.map(qseq,MD=True,cs=True))

    if best_n > 1:
        return itr

    if not itr:
        return

    for hit in itr:
        #TODO if hit.r_st==intron_len <- look for additional alignments?
        #or it is not needed because the longest alignment is the "best" one?

        if hit.strand != 1:
            continue

        if not hit.is_primary:
            continue

        else:
            return hit #returns the first primary hit

def calc_alignment_score(hit,scoring):
    if hit is None:
        return None

    # TODO verify alignment_score calc
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

    alignment_score = matches*scoring[0] - (mismatches)*scoring[1] - gap_penalty

    return alignment_score

def get_alignment(chrom, jstart, jend, overhang, pysam_fa, max_length, scoring,  k, w, m):


    intron_len = jend - jstart

    rstart = max(jstart - overhang, 0)
    rend = min(jstart + overhang, max_length)
    qstart = max(jend - overhang, 0)
    qend = min(jend + overhang, max_length)

    rseq = get_seq(chrom, rstart, rend, pysam_fa)
    qseq = get_seq(chrom, qstart, qend, pysam_fa )
    hit = align_seq_pair(rseq, qseq, scoring,k,w,m)

    if hit:
        #check if the alignment is in the overlap region of short introns
        if overhang * 2 >= intron_len:

            if hit.r_st == intron_len:
                return None

            if overhang > intron_len:
                if hit.r_st - hit.q_st == intron_len:
                    return None

            if hit.r_st >= overhang and hit.q_en <= overhang:
                return None

    return hit


# def get_self_alignment(chrom, start, end, max_length, scoring, k, w, m, pysam_fa,  overhang):
#     rseq,qseq,hit = get_alignment(chrom, start, end, overhang, pysam_fa, max_length, scoring,  k, w, m)
#     return rseq,qseq,hit


def get_flanking_subsequences(introns,chrom_sizes,overhang,ref_fa):
    tmp_regions = tempfile.NamedTemporaryFile(mode='a',dir=os.getcwd(),delete=False)
    tmp_fa = tempfile.NamedTemporaryFile(dir=os.getcwd(),delete=False)

    seen = set()
    for key in list(introns.keys()):
        chrom = key[0]
        jstart = key[1]
        jend = key[2]
        max_length = chrom_sizes[chrom]
        rstart = max(jstart - overhang + 1, 1) #1 based
        rend = min(jstart + overhang, max_length)
        qstart = max(jend - overhang + 1, 1) #1 based
        qend = min(jend + overhang, max_length)
        r1 = f'{chrom}:{rstart}-{rend}'
        r2 = f'{chrom}:{qstart}-{qend}'
        introns[key]['jstart'] = r1
        introns[key]['jend'] = r2

        if rstart> max_length or qstart > max_length:
            #remove key from introns dict
            del introns[key]
            print(f'Warning: intron {key} from the GTF file is out of range in FASTA file. Skipping...')
            continue

        for r in [r1,r2]:
            if r not in seen:
                t=tmp_regions.write(f'{r}\n')
                seen.add(r)


    tmp_regions.close()
    tmp_fa.close()

    cmd1 = f"samtools faidx {ref_fa} -r {tmp_regions.name} -o {tmp_fa.name}"
    p1 = subprocess.run(shlex.split(cmd1), check=True)


    seqs = {}
    with open(tmp_fa.name,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                seq_name = line[1:]
                if seq_name not in seqs:
                    seqs[seq_name] = ''
                continue
            sequence = line
            seqs[seq_name]=seqs[seq_name] + sequence

    os.unlink(tmp_regions.name)
    os.unlink(tmp_fa.name)

    return seqs
