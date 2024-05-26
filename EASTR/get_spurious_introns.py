import collections
import os
import shlex
import subprocess
import tempfile
import time

import pysam

from EASTR import extract_junctions
from EASTR import alignment_utils
from EASTR import utils


def get_self_aligned_introns(introns, seqs, overhang, k, w, m, scoring):
    self_introns = {}
    for key, value in introns.items():
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


def linear_distance(string1, string2):
    if len(string1)!=len(string2):
        raise ValueError("strings must be of equal length")
    distance = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            distance += 1
    return distance

def get_middle_seq(len_,seq):
    if len(seq) > len_:
        middle_start = (len(seq) - len_) // 2
        middle_end = middle_start + len_
        seq = seq[middle_start:middle_end]
    return seq


def bowtie2_align_self_introns_to_ref (introns_to_align, seqs, bt2_index, overhang, p=1, len_=15, bt2_k=10):
    bt2_k = bt2_k + 1
    tmp_sam = tempfile.NamedTemporaryFile(dir=os.getcwd(),delete=False)
    tmp_fa = tempfile.NamedTemporaryFile(mode='a',dir=os.getcwd(),delete=False)
    for key,value in introns_to_align.items():
        read_name = ','.join([str(i) for i in key])
        hit = value['hit']
        rseq = seqs[value['jstart']]
        qseq = seqs[value['jend']]
        seq1 = get_middle_seq(len_*2, rseq[hit.r_st:hit.r_en])
        seq2 = get_middle_seq(len_*2, qseq[hit.q_st:hit.q_en])
        seqh = rseq[overhang-len_:overhang] + qseq[overhang:overhang+len_]

        x = tmp_fa.write(f'>{read_name},seq1\n{seq1}\n' + \
                         f'>{read_name},seq2\n{seq2}\n' + \
                         f'>{read_name},seqh\n{seqh}\n')
    tmp_fa.close()
    tmp_sam.close()

    #TODO -R {bt2_k} -N 2?
    cmd = f"bowtie2 -p {p} --end-to-end -k {bt2_k} -D 20 -R 5 -L 20 -N 1 -i S,1,0.50 -x {bt2_index} -f {tmp_fa.name} -S {tmp_sam.name}"
    subprocess.run(shlex.split(cmd),stderr=subprocess.DEVNULL, check=True)
    samfile = pysam.AlignmentFile(tmp_sam.name,'r')

    d = collections.defaultdict(list)
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

def is_two_anchor_alignment(hit,overhang,anchor):
    c1 = (hit.r_en < overhang + anchor - 1)
    c2 = (hit.r_st > overhang - anchor)
    c3 = (hit.q_en < overhang + anchor - 1)
    c4 = (hit.q_st > overhang - anchor)
    c5 = (abs(hit.r_st - hit.q_st) > anchor * 2)
    if (c1 or c2 or c3 or c4 or c5):
        return False
    else:
        return True


def is_spurious_alignment(key, value, seqs, overhang, min_duplicate_exon_length, bt2_k=10, anchor=7):
    is_two_anchor = is_two_anchor_alignment(value['hit'],overhang,anchor)
    hit = value['hit']

    if is_two_anchor:
        #check unique alignment
        if (value['seq1'] == 1) or (value['seq2'] == 1):
            if value['seqh'] == 0:
                # print("unique alignment")
                # print(key)
                return False


    else:
        #if duplicated exon:
        if hit.q_st - hit.r_st >= min_duplicate_exon_length:
            if value['seqh'] == 0:
                # print("duplicated exon")
                # print(key)
                return False


        if (value['seq1'] < bt2_k ) or (value['seq2'] < bt2_k):
            if value['seqh'] == 0:
                # print("partial alignment - no hybrid sequence found")
                # print(key)
                return False

            rseq = seqs[value['jstart']]
            qseq = seqs[value['jend']]
            e5 = rseq[overhang - anchor:overhang]
            i5 = rseq[overhang:overhang + anchor]
            i3 = qseq[overhang - anchor:overhang]
            e3 = qseq[overhang:overhang + anchor]
            distance_e5i3 = linear_distance(e5, i3)
            distance_i5e3 = linear_distance(i5, e3)
            if distance_e5i3 + distance_i5e3 > 2:
                # print("partial alignment - no overhang")
                # print(key)
                return False
            else:
                # print("partial alignment - has overhang")
                # print(key)
                return True

    return True


def get_spurious_introns(self_introns, seqs, bt2_index, overhang, min_duplicate_exon_length,
                          anchor=7, min_junc_score=1, p=1, is_bam=True, bt2_k=10):
    introns_to_align = {}
    spurious = {}
    for k,v in self_introns.items():

        if is_bam and v['score'] <= min_junc_score:
            spurious[k] = v
            continue

        introns_to_align[k] = v

    bw2_alignments = bowtie2_align_self_introns_to_ref(introns_to_align, seqs, bt2_index, overhang, p=p, len_=15, bt2_k=bt2_k)

    for k, v in introns_to_align.items():
        if is_spurious_alignment(k, v, seqs, overhang, min_duplicate_exon_length, anchor=anchor):
            spurious[k] = v

    return spurious


def get_spurious_junctions(scoring, k, w, m, overhang, min_duplicate_exon_length, bt2_index, bt2_k, ref_fa, p, anchor,
                            min_junc_score, bam_list, gtf_path,
                                    bed_path, trusted_bed, out_original_junctions, verbose):
    chrom_sizes = utils.get_chroms_list_from_fasta(ref_fa)

    if bam_list:
        is_bam = True

        #make tmp files if not given out_original_junctions
        remove_tmp=False
        utils.make_dir('tmp')
        if out_original_junctions is None:
            remove_tmp = True
            out_original_junctions = [tempfile.NamedTemporaryFile(delete=False, dir='tmp', suffix='.bed') for i in range(len(bam_list))]
            #close all the temporary files in out_original_junctions:
            for f in out_original_junctions:
                f.close()
            out_original_junctions = [f.name for f in out_original_junctions]


        if verbose:
            start_extr = time.time()
            print('extracting junctions from bam files...')

        introns = extract_junctions.junction_extractor_multi_bam(bam_list,out_original_junctions,p)

        if verbose:
            end_extr = time.time()
            print('extracting junctions took {} seconds'.format(end_extr-start_extr))

        #delete all the temporary files in out_original_junctions:
        if remove_tmp:
            for f in out_original_junctions:
                os.remove(f)

    elif gtf_path:
        is_bam = False
        introns = extract_junctions.extract_splice_sites_gtf(gtf_path)

    elif bed_path:
        is_bam = False
        if len(bed_path) == 1:
            introns = extract_junctions.get_junctions_multi_bed(bed_path,p)


    else:
        raise ValueError('No input file given')

    if trusted_bed:
        trusted_introns = extract_junctions.get_junctions_from_bed(trusted_bed)
        introns = {k:v for k,v in introns.items() if k not in trusted_introns}

    if verbose:
        print('Getting spurious junctions...')
        start_spur = time.time()

    seqs = alignment_utils.get_flanking_subsequences(introns, chrom_sizes, overhang, ref_fa)
    self_introns = get_self_aligned_introns(introns, seqs, overhang, k, w, m, scoring)
    spurious_dict = get_spurious_introns(self_introns, seqs, bt2_index, overhang, min_duplicate_exon_length, anchor=anchor,
                                     min_junc_score=min_junc_score, p=p, is_bam=is_bam, bt2_k=bt2_k)
    spurious_dict  = dict(sorted(spurious_dict.items(), key=lambda x: (x[0][1], x[0][1], x[0][2], x[0][3])))

    if verbose:
        end_spur = time.time()
        print('Getting spurious junctions took {} seconds'.format(end_spur - start_spur))

    return spurious_dict
