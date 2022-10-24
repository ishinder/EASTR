
from typing import Dict
import pandas as pd
import pysam
import mappy as mp
import collections
import os
import subprocess
import shlex
from io import StringIO
import re
from EASTR import utils, filter_bam
from EASTR.alignment_utils import *
import numpy as np

class Alignment:
    def __init__(self,alignment):
        self.query_name = alignment.query_name
        self.reference_name = alignment.reference_name
        self.reference_start = alignment.reference_start
        self.NH = alignment.get_tag("NH")
        self.is_proper_pair = alignment.is_proper_pair
        self.is_read1 = alignment.is_read1
        self.next_reference_start = alignment.next_reference_start


def get_introns_from_bam(samfile):
    introns=collections.OrderedDict() 
    for alignment in samfile.fetch(until_eof=True, multiple_iterators=True):
        currentloc = alignment.pos
        for i,cigarop in enumerate(alignment.cigar):
            if (cigarop[0]==4): #substitution or insertion in query
                continue
            if (cigarop[0]==1):
                continue
            if(cigarop[0]==3):
                key = (alignment.reference_name,currentloc,currentloc+cigarop[1])
                #TODO check that previous + next have cigarop==0
                o5 = alignment.cigar[i-1][1]           
                o3 = alignment.cigar[i+1][1]
                a = Alignment(alignment)
                if key not in introns:
                    introns[key]=[(a,o5,o3)]
                else:
                    introns[key].append((a,o5,o3))  
                    #TODO - first o5 is longest? last o3 is longest?
                    #TODO - check if any key has duplicate alignments??? len(set)==len(list)          
            currentloc=currentloc+cigarop[1]
    return introns
    

def find_spurious_alignments(introns, samfile, ref_fa, chrom_sizes, read_length, scoring, k, w, min_chain_score, filter_bam=True):
    spurious_introns = {}
    spurious_alignments = set()
    NH = collections.defaultdict((int))
    seen_alignments = set()
    
    for key, values in introns.items():
        max_length = chrom_sizes[key[0]]
        o5 = max(values, key=lambda x: x[1])[1]
        o3 = max(values, key=lambda x: x[2])[2]
        hits = get_alignment(key[0], key[1], key[2], o5, o3,
                             ref_fa, max_length, read_length, scoring,  k, w,min_chain_score)
        
        score = max(hits, key=lambda x: x[1])[1]
        if score > -scoring[1]*read_length:    
            spurious_introns[key] = score 
            #TODO: any additional info to include -reads, alignments, etc?
            
            if not filter_bam: #do not need to find spurious alignments if a filtered bam is not output. 
                continue 

            alignments = [x[0] for x in values]
            for alignment in alignments:
                akey = (alignment.query_name, alignment.reference_name, alignment.reference_start)
                if akey not in seen_alignments:
                    seen_alignments.add(akey)
                    spurious_alignments.add(akey)
                    NH[(alignment.query_name,alignment.is_read1)] += 1

                    if alignment.is_proper_pair:
                        mkey = (alignment.query_name, alignment.reference_name, alignment.next_reference_start)
                        seen_alignments.add(mkey)
                        spurious_alignments.add(mkey)
                        NH[(alignment.query_name, not alignment.is_read1)] += 1

    return spurious_alignments, spurious_introns, NH

def write_filtered_bam(outbam, samfile, spurious_alignments, NH):
    removed_reads = set()
    outf = pysam.AlignmentFile(outbam, "wb", template=samfile)

    for alignment in samfile.fetch(until_eof=True, multiple_iterators=True):
        if alignment.is_unmapped: #TODO - decide if unmapped alignments should be included in output bam file (to convert back to fastq)?
            continue 

        key = (alignment.query_name, alignment.reference_name, alignment.reference_start)
        if key in spurious_alignments:
            if alignment.get_tag('NH') == NH[(alignment.query_name, alignment.is_read1)]:
                removed_reads.add((alignment.query_name, alignment.is_read1))
            continue
            #TODO: add custom tag instead of removing?: alignment.tags = alignment.tags + [('XR',1)]
            #alignment.is_mapped = False
            # alignment.mapq = 0
            # alignment.cigar = []
            # alignment.cigarstring=None
            
            
        
        if (alignment.query_name, alignment.is_read1) in NH:
            new_NH = alignment.get_tag("NH") - NH[(alignment.query_name,alignment.is_read1)]
            # print(alignment.qname,alignment.get_tag("NH"),new_NH)
            alignment.set_tag("NH", new_NH)
            if new_NH == 0: #TODO remove? (this should not happen)
                os.remove(outbam)
                raise Exception("NH tag cannot be zero")
                # print("removed read: ",alignment.qname,alignment.is_read1)
                continue

            if new_NH < 0:  #TODO remove? (this should not happen)
                os.remove(outbam)
                raise Exception("NH tag cannot be negative")


        w=outf.write(alignment) 
    outf.close()
    return removed_reads



def filter_alignments_from_bam(ref_fa, bam, scoring, read_length, k, w, m, outbam=None):

    #check if file is cram or bam:
    if bam.split('.')[-1]=="bam":
        samfile = pysam.AlignmentFile(bam, "rb") # type: ignore
    
    else :
        samfile = pysam.AlignmentFile(bam, "rc") 

    removed_reads = set()
    chrom_sizes = utils.get_chroms_list_from_bam(bam)
    
    introns = get_introns_from_bam(samfile)
    spurious_alignments, spurious_introns, NH = find_spurious_alignments(introns, samfile, ref_fa,
                                                                        chrom_sizes, read_length, scoring, k, w, m)


    #write new bam
    if outbam is not None:
        removed_reads = write_filtered_bam(outbam, samfile, spurious_alignments, NH)
    
    return spurious_introns, removed_reads


    

if __name__ == '__main__':
    import time
    scoring=[3,3,4,2,24,1,1]
    
    k = 15
    w = 10
    m = 15
    
    start = time.time()
    bam = "tests/data/ERR188044_chrX.bam"
    read_length = utils.get_read_length_from_bam(bam)
    ref_fa = "tests/data/chrX.fa"
    outbam = "tests/output/ERR188044_chrX_filtered.bam"

    spur_introns,removed_reads = filter_alignments_from_bam(ref_fa, bam, scoring, read_length, k, w,m,outbam=outbam)
    end = time.time()
    print(f"took {end-start} seconds")
    
    # bam="/ccb/salz8-2/shinder/projects/EASTR_tests/chess_brain/BAM/original/R2824_C4KHUACXX.bam"
    # outbam = "/ccb/salz8-2/shinder/projects/EASTR_tests/chess_brain/BAM/filtered/R2824_C4KHUACXX_A3B3k15m15w10v5_filtered.bam"
    # read_length=utils.get_read_length_from_bam(bam)
    # ref_fa= "/ccb/salz8-2/chess-brain/ref/hg38mod_noPARs.fa"
    # spur_introns,removed_reads = filter_alignments_from_bam(ref_fa, bam, scoring, read_length, k, w,m,outbam=outbam)

    # bam="/ccb/salz8-2/shinder/projects/EASTR_tests/GTEx2/BAM/original/GTEX-11H98-0011-R11b-SM-5NQ6U.sorted.bam"
    # read_length = utils.get_read_length_from_bam(bam)
    # ref_fa= "/ccb/salz8-2/chess-brain/ref/hg38mod_noPARs.fa"
    # outbam = "/ccb/salz8-2/shinder/projects/EASTR_tests/GTEx2/BAM/filtered/GTEX-11H98-0011-R11b-SM-5NQ6U_A3B3k15m15w10v5_filtered.bam"

    # start = time.time()
    # bam = "/ccb/salz8-2/shinder/projects/Geuvadis/BAM/ERR188025.bam"
    # outbam = "/ccb/salz8-2/shinder/projects/Geuvadis/EASTR2/BAM/ERR188025_filtered.bam"
    # out_introns = "/ccb/salz8-2/shinder/projects/Geuvadis/EASTR2/BED/ERR188025.bed"
    # read_length = utils.get_read_length_from_bam(bam)
    # outdir="/ccb/salz8-2/shinder/projects/Geuvadis/EASTR/BAM"
    # ref_fa= "/ccb/salz8-2/chess-brain/ref/hg38mod_noPARs.fa"
    # spurious_introns, removed_reads = filter_alignments_from_bam(ref_fa, bam, scoring, read_length, k, w, m, outbam=outbam)
    # end = time.time()
    # print(f"took {(end-start)/60} minutes")

    # #chroms = [x['SN'] for x in samfile.header['SQ']]
    
    # start = time.time()
    # bam = "/ccb/salz8-2/shinder/projects/repeat_spliced_alignments/chess3_tiebrush_bam/all_tissues.tb.bam"
    # introns = get_spurious_junctions_from_bam(bam, scoring, ref_fa, read_length,k,w)
    # end = time.time()
    # print(f"took {end-start} seconds")
    
    # introns[~introns['score'].isna()]