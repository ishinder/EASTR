import csv
from posixpath import basename
import re
from EASTR import extract_junctions
from EASTR import alignment_utils
from EASTR.utils import get_chroms_list_from_fasta
# import mappy as mp
import pandas as pd
import subprocess
import pysam
import shlex
from collections import defaultdict
import os
import multiprocessing
import tempfile




# def remove_low_score_junctions(introns,trusted_junctions={},n=0):
#     if n == 0:
#         return introns
#     low_score_junctions = set()
#     for key in introns:
#         if introns[key]['score'] <= n:
#             if key in trusted_junctions:
#                 continue
#             low_score_junctions.add(key)

#     filtered_introns = introns
#     for key in introns:
#         if key in low_score_junctions:
#             filtered_introns.pop(key)

#     return filtered_introns



def get_self_aligned_introns(introns, seqs, overhang, k, w, m, scoring):
    self_introns = {}
    for key,value in introns.items():
        rseq = seqs[introns[key]['jstart']]
        qseq = seqs[introns[key]['jend']]
        hit = alignment_utils.align_seq_pair(rseq, qseq, scoring,k,w,m)
        if hit:
            intron_len = key[2] - key[1]
            if overhang * 2 >= intron_len:
                if hit.r_st == intron_len:
                    hit = None
                elif overhang > intron_len:
                    if hit.r_st - hit.q_st == intron_len:
                        hit = None
                elif hit.r_st >= overhang and hit.q_en <= overhang:
                    hit = None
        if hit:
            self_introns[key] = value
            self_introns[key]['hit'] = hit
    return self_introns

def get_introns_multi_bam_regtools(bam_list,p):
    pool = multiprocessing.Pool(p)
    results = pool.map(extract_junctions.extract_splice_sites_bam_regtools,bam_list)
    pool.close
    dd = defaultdict(dict)
    for i,d in enumerate(results):
        sample_id = re.split(r'.[bcr]+am',basename(bam_list[i]))[0]
        for key, value in d.items():
            
            if key not in dd:
                dd[key]['samples']= set([sample_id])
                dd[key]['score'] = value[0]
                dd[key]['score_indiv_max'] = value[0]
                dd[key]['o5'] = value[1]
                dd[key]['o3'] = value[2]
            else:
                dd[key]['samples'].add(sample_id)
                dd[key]['score'] = dd[key]['score'] + value[0]
                dd[key]['score_indiv_max'] = max(dd[key]['score_indiv_max'],value[0])
                dd[key]['o5'] = max(dd[key]['o5'],value[1])
                dd[key]['o3'] = max(dd[key]['o3'],value[2])
    return dd


# def get_introns_multi_bam_custom(bam_list,p):
#     pool = multiprocessing.Pool(p)
#     results = pool.map(extract_junctions.extract_splice_sites_bam_custom,bam_list)
#     pool.close
#     dd = defaultdict(dict)
#     for i,d in enumerate(results):
#         for key, value in d.items():
#             if key not in dd:
#                 dd[key]['bam_path'] = [{bam_list[i]:}]
#                 dd[key]['score'] = 
#                 dd[key]['o5'] = 
#                 dd[key]['o3'] = 
#             else:
#                 dd[key]['bam_path'].append(bam_list[i])
#                 dd[key]['score'] += value[0]
#                 dd[key]['o5'] = max(dd[key]['o5'],value[1])
#                 dd[key]['o3'] = max(dd[key]['o3'],value[2])
#     return dd



def linear_distance(string1, string2):
    if len(string1)!=len(string2):
        raise Exception("strings must be of equal length")
    distance = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            distance += 1
    return distance


def bowtie2_align_self_introns_to_ref (introns_to_align, seqs, bt2_index, overhang, p=1, hseq_len=25, bt2_k=10):
    bt2_k = bt2_k + 1
    tmp_sam = tempfile.NamedTemporaryFile(dir=os.getcwd(),delete=False)
    tmp_fa = tempfile.NamedTemporaryFile(mode='a',dir=os.getcwd(),delete=False)
    for key,value in introns_to_align.items():
        read_name = ','.join([str(i) for i in key])
        hit = value['hit']
        rseq = seqs[value['jstart']]
        qseq = seqs[value['jend']]
        seq1 = rseq[hit.r_st:hit.r_en]
        seq2 = qseq[hit.q_st:hit.q_en]
        seqh = rseq[-hseq_len:] + qseq[overhang:overhang+hseq_len]

        x = tmp_fa.write(f'>{read_name},seq1\n{seq1}\n' + \
                         f'>{read_name},seq2\n{seq2}\n' + \
                         f'>{read_name},seqh\n{seqh}\n')
    tmp_fa.close()
    tmp_sam.close()

    #TODO -R {bt2_k} -N 2? 
    cmd = f"bowtie2 -p {p} --end-to-end -k {bt2_k} -D 20 -R 5 -L 20 -N 1 -i S,1,0.50 -x {bt2_index} -f {tmp_fa.name} -S {tmp_sam.name}"
    subprocess.run(shlex.split(cmd),stderr=subprocess.DEVNULL)
    samfile = pysam.AlignmentFile(tmp_sam.name,'r')
    
    d = defaultdict(list)
    for alignment in samfile.fetch(until_eof=True):
        qname= alignment.qname.split(',')
        qname[1] = int(qname[1])
        qname[2] = int(qname[2])
        qname = tuple(qname)
        d[qname].append(alignment)

        if qname[4] == 'seqh':
            if alignment.is_unmapped:
                introns_to_align[qname[0:4]][qname[4]] = 0
                continue

        if qname[4] not in introns_to_align[qname[0:4]]:
            introns_to_align[qname[0:4]][qname[4]] = 1
        
        else:
            introns_to_align[qname[0:4]][qname[4]] += 1

    os.unlink(tmp_sam.name)
    os.unlink(tmp_fa.name)

    return d #TODO return?


def is_spurious_alignment(key, value, seqs, overhang, bt2_k=10, anchor=7):
    hit = value['hit']

    #check unique junction
    if hit.q_st - hit.r_st > overhang:
        if value['seqh'] == 0:
            return False

    #check unique alignment
    if (value['seq1'] == 1) and (value['seq2'] == 1):
        if value['seqh'] == 0:
            return False

    #check other conditions
    c1 = (hit.r_en <= overhang + anchor)
    c2 = (hit.r_st >= overhang - anchor)
    c3 = (hit.q_en <= overhang + anchor)
    c4 = (hit.q_st >= overhang - anchor)
    c5 = (abs(hit.r_st - hit.q_st) > anchor * 2)

    if (c1 or c2 or c3 or c4 or c5):
        if (value['seq1'] < bt2_k ) or (value['seq2'] < bt2_k):
            rseq = seqs[value['jstart']]
            qseq = seqs[value['jend']]
            e5 = rseq[overhang - anchor:overhang]
            i5 = rseq[overhang:overhang + anchor]
            i3 = qseq[overhang - anchor:overhang]
            e3 = qseq[overhang:overhang + anchor]
            distance_e5i3 = linear_distance(e5, i3)
            distance_i5e3 = linear_distance(i5, e3)
            if distance_e5i3 + distance_i5e3 > 2:
                return False

    return True



def get_spurious_introns(self_introns, seqs, bt2_index, overhang, anchor=7, min_junc_score_global=1, 
                        min_junc_score_local=1, p=1, trusted_introns={}, is_bam=True, bt2_k=10):
    introns_to_align = {}
    spurious = {}
    for k,v in self_introns.items():
        if k in trusted_introns:
            continue

        if is_bam and (v['score'] <= min_junc_score_global and v['score_indiv_max'] <= min_junc_score_local):  
            spurious[k] = v
            continue
        
        introns_to_align[k] = v
 
    bw2_alignments = bowtie2_align_self_introns_to_ref(introns_to_align, seqs, bt2_index, overhang, p=p, hseq_len=25, bt2_k=bt2_k)

    for k, v in introns_to_align.items():
        if is_spurious_alignment(k, v, seqs, overhang, anchor=anchor):
            spurious[k] = v
    
    return spurious

# def spurious_junctions_in_bam(bam_list, scoring, k, w, m, overhang, bt2_index, ref_fa, p, anchor, min_junc_score, trusted_introns={}):
#     introns = get_introns_multi_bam(bam_list,p)
#     chrom_sizes = get_chroms_list_from_fasta(ref_fa)
#     seqs = alignment_utils.get_flanking_subsequences(introns,chrom_sizes,overhang,ref_fa)
#     self_introns = get_self_aligned_introns(introns, seqs, overhang, k, w, m, scoring)
#     spurious = get_spurious_introns(self_introns, seqs, bt2_index, overhang, anchor=7, min_junc_score=min_junc_score, p=p, trusted_introns={})
#     return spurious


def spurios_junctions_in_gtf(gtf_path, scoring, k, w, m, overhang, bt2_index, ref_fa, p, anchor, min_junc_score):
    introns = extract_junctions.extract_splice_sites_gtf(gtf_path)
    chrom_sizes = get_chroms_list_from_fasta(ref_fa)
    seqs = alignment_utils.get_flanking_subsequences(introns,chrom_sizes,overhang,ref_fa)
    self_introns = get_self_aligned_introns(introns, seqs, overhang, k, w, m, scoring)
    spurious = get_spurious_introns(self_introns, seqs, bt2_index, overhang, anchor=7, min_junc_score_global=min_junc_score, p=p, is_bam=False)
    return spurious

def get_spurious_junctions(scoring, k, w, m, overhang, bt2_index, bt2_k, ref_fa, p, anchor, 
                            min_junc_score_global, min_junc_score_local, **kargs):
    chrom_sizes = get_chroms_list_from_fasta(ref_fa)

    if kargs['trusted_introns'] is not None:
        trusted_introns = kargs['trusted_introns']
    else:
        trusted_introns = {}

    if kargs['bam_list'] is not None:
        bam_list = kargs['bam_list']
        is_bam = True
        introns = get_introns_multi_bam_regtools(bam_list,p)

    else: 
        gtf_path = kargs['gtf_path']
        is_bam = False
        introns = extract_junctions.extract_splice_sites_gtf(gtf_path)

    seqs = alignment_utils.get_flanking_subsequences(introns, chrom_sizes, overhang, ref_fa)
    self_introns = get_self_aligned_introns(introns, seqs, overhang, k, w, m, scoring)
    spurious = get_spurious_introns(self_introns, seqs, bt2_index, overhang, anchor=anchor, min_junc_score_global=min_junc_score_global,
                                     min_junc_score_local=min_junc_score_local, p=p, trusted_introns=trusted_introns, is_bam=is_bam,
                                     bt2_k=bt2_k)

    return spurious

def spurious_junctions_to_bed(spurious, fileout, is_bam = True, write=True, **kargs):
    df = pd.Series(spurious).reset_index()
    l = len(str(len(spurious)))
    df['name'] = df.index + 1
    df['name'] = df['name'].apply(lambda x: f'{{0:0>{l}}}'.format(x))
    df['name'] = 'JUNC' + df['name']
   
    if is_bam:
        #for bam, the score is the count of reads across all samples supporting a given junction
        df['score'] = pd.DataFrame(df[0].tolist())['score']
        df['samples'] = pd.DataFrame(df[0].tolist())['samples']
        df['samples'] = [','.join(map(str, l)) for l in df['samples']]
        df = df.iloc[:,list(range(3))+ [-2, -3, 3,-1]]
    else:
        scoring = kargs['scoring']
        #for gtf, the score is the alignment score (per scoring matrix) of a given junction
        df['transcripts'] = pd.DataFrame(df[0].tolist())['transcripts']
        df['transcripts'] = [','.join(map(str, l)) for l in df['transcripts']]
        
        for i, row in df.iterrows():
            key = tuple(row[0:4])
            hit = spurious[key]['hit']
            score = alignment_utils.calc_alignment_score(hit, scoring)
            df.at[i,'score'] = score
        df['score'] = df['score'].astype(int)
        df = df.iloc[:,list(range(3)) + [-3, -1, 3, -2]]

    if write:
        df.to_csv(fileout, header=False, index=False, sep='\t', quoting=csv.QUOTE_NONE)

    return df


if __name__ == '__main__':
    import time
    import glob

    scoring=[3,4,12,2,32,1,1]
    k = 3
    w = 2
    m = 25
    overhang = 50
    bt2_k = 10
    min_junc_score_global = 1
    min_junc_score_local = 1
    anchor = 7
    p = 18
    ref_fa = "/ccb/salz7-data/genomes/hg38/hg38p13.fa"
    bt2_index = '/ccb/salz8-2/shinder/projects/EASTR/tests2/other/hg38mod_noPARs_index/hg38mod_noPARs_index'
    trusted_introns = None
    is_bam = True
    gtf_path = None

    #arabidopsis
    ref_fa="/ccb/salz8-3/shinder/projects/EASTR_tests/arabidopsis/annotation/GCA_000001735.2_TAIR10.1_genomic.fna"
    bt2_index="/ccb/salz8-3/shinder/projects/EASTR_tests/arabidopsis/annotation/bt2_GCA_000001735.2_TAIR10.1_genomic/GCA_000001735.2_TAIR10.1_genomic"
    bam_list=glob.glob("/ccb/salz8-3/shinder/projects/EASTR_tests/arabidopsis/BAM/original/*bam")
    gtf_path= "/ccb/salz8-3/shinder/projects/EASTR_tests/arabidopsis/annotation/GCA_000001735.2_TAIR10.1_genomic.gtf"
    basedir='/ccb/salz8-3/shinder/projects/EASTR_tests/arabidopsis'
    

    #maize
    ref_fa="/ccb/salz2/shinder/projects/EASTR_tests/maize/ref/B73/Zm-B73-REFERENCE-NAM-5.0.fa"
    bam_list=glob.glob("/ccb/salz2/shinder/projects/EASTR_tests/maize/BAM/original/*d0.1*bam")
    bt2_index="/ccb/salz2/shinder/projects/EASTR_tests/maize/ref/B73/bowtie2/B73-NAM-5.0"
    spurious=get_spurious_junctions(scoring, k, w, m, overhang, bt2_index, bt2_k, ref_fa, p, anchor, 
                            min_junc_score_global, min_junc_score_local, bam_list=bam_list, trusted_introns=None)
    out_dir="/ccb/salz2/shinder/projects/EASTR_tests/maize/output"
    out_introns=f"{out_dir}/spurious_introns_in_2_d.01_samples.bed"
    suffix=".filtered"
    is_bam=True


    # bam_list = glob.glob('/ccb/salz8-2/shinder/projects/EASTR_tests/chess_brain/BAM/original/*.bam')

    bam_list = glob.glob('/ccb/salz8-2/shinder/projects/EASTR/tests/data/chrX_data/hisat2/*.bam')
    bam_list = "/ccb/salz8-2/shinder/projects/EASTR/tests/data/chrX_data/bamlist.txt"
    gtf_path = '/ccb/salz8-2/shinder/projects/EASTR/EASTR/tests/data/chrX_data/genes/chrX.gtf'

    out_introns = "/ccb/salz8-2/shinder/projects/EASTR/EASTR/tests/output/chrX_data/spurious_introns.bed"
    out_dir = "/ccb/salz8-2/shinder/projects/EASTR/tests/output/chrX_data/eastr"


    
#     spurious = spurious_junctions_in_bam(bam_list, scoring, k, w, m, overhang, bt2_index, ref_fa, p, anchor, min_junc_score_global)

#     # gtf_path = "/ccb/salz8-2/shinder/projects/EASTR/tests/data/chess3.0.gtf"
    
#     fileout = "/ccb/salz8-3/shinder/misc/for_KH/spurious_introns_from_n15_riboZ_CHESS_brain.bed"

#     ref_fa_gtf = '' #ref
    


        
#     gtf_path = "/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf"

    
    
#     start = time.time()
    
#     end = time.time()

#     bw2_alignments = bowtie2_align_self_introns_to_ref(self_introns, seqs, bt2_index, overhang, p=1, hseq_len=25)
#     print(f'{end-start}/60 seconds')