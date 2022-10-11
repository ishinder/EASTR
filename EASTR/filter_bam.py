
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


def get_mate(alignment,index):
    if alignment.mate_is_unmapped:
        return None
    mate_start = alignment.next_reference_start
    mate_chrom = alignment.next_reference_name
    tlen = alignment.tlen
    for mate in index.find(alignment.qname):
        if ((mate_start== mate.reference_start) & (mate_chrom == mate.reference_name) & (tlen == -mate.tlen)):
            return mate



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
                key = (samfile.getrname(alignment.tid),currentloc,currentloc+cigarop[1])
                o5 = alignment.cigar[i-1][1]           
                o3 = alignment.cigar[i+1][1]
                if key not in introns:
                    introns[key]=[(alignment,o5,o3)]
                else:
                    introns[key].append((alignment,o5,o3))  
                    #TODO - first o5 is longest? last o3 is longest?
                    #TODO - check if any key has duplicate alignments??? len(set)==len(list)          
            currentloc=currentloc+cigarop[1]
    return introns
    
def index_samfile_by_reads(samfile):
    index = pysam.IndexedReads(samfile)
    index.build()
    return index

def find_spurious_alignments(introns, index, ref_fa, chrom_sizes, read_length, scoring, k, w, min_chain_score):
    spurious_introns = {}
    spurious_alignments = set()
    spurious_mates = set()
    NH = collections.defaultdict(lambda:0)
    seen_alignments = set()
    
    for key, values in introns.items():
        max_length = chrom_sizes[key[0]]
        o5 = max(values, key=lambda x: x[1])[1]
        o3 = max(values, key=lambda x: x[2])[2]
        hits = get_alignment(key[0], key[1], key[2], o5, o3,
                             ref_fa, max_length, read_length, scoring,  k, w,min_chain_score)
        
        score = max(hits, key=lambda x: x[1])[1]
        if score > -scoring[1]*read_length:    
            spurious_introns[key]=score 
            #TODO: any additional info to include -reads, alignments, etc?
            
            alignments = [x[0] for x in values]
            for alignment in alignments:
                if alignment not in seen_alignments:
                    NH[(alignment.qname,alignment.is_read1)] += 1
                    mate = get_mate(alignment, index)
                    seen_alignments.add(alignment)
                    spurious_alignments.add(alignment)
                    spurious_mates.add(mate)

    return spurious_alignments, spurious_mates, spurious_introns, NH

def write_filtered_bam(outbam, samfile, spurious_alignments, spurious_mates, NH, index):
    removed_reads = set()
    outf = pysam.AlignmentFile(outbam, "wb", template=samfile)

    for alignment in samfile.fetch(until_eof=True, multiple_iterators=True):
        if alignment.is_unmapped:
            continue 
        if alignment in spurious_mates: #this is a mate of a spurious alignment
            if alignment.is_proper_pair: #throwing out the second pair improves precision + sensitivity
                removed_reads.add((alignment.qname,alignment.is_read1))
                continue

            alignment.mate_is_mapped = False
            alignment.tlen = 0
            alignment.next_reference_start = alignment.reference_start
            alignment.next_reference_name = alignment.reference_name
            alignment.next_reference_id  = alignment.reference_id
            #TODO check if any other tags/fields need to be updated

        if alignment in spurious_alignments:
            continue
            #TODO: add custom tag instead?: alignment.tags = alignment.tags + [('XR',1)]
            #alignment.is_mapped = False
            # if alignment.is_read2:
            #     alignment.flag=133
            # else:
            #     alignment.flag=69
            # alignment.mapq = 0
            # alignment.cigar = []
            # alignment.cigarstring=None
            
            
        
        if (alignment.qname, alignment.is_read1) in NH:
            new_NH = alignment.get_tag("NH") - NH[(alignment.qname,alignment.is_read1)]
            # print(alignment.qname,alignment.get_tag("NH"),new_NH)
            alignment.set_tag("NH", new_NH)
            if new_NH == 0:
                removed_reads.add((alignment.qname,alignment.is_read1))
                # print("removed read: ",alignment.qname,alignment.is_read1)
                continue

            if new_NH < 0:
                os.remove(outbam)
                raise Exception("NH tag cannot be negative")


        w=outf.write(alignment) 
    outf.close()
    return removed_reads



def filter_alignments_from_bam(ref_fa, bam, scoring, read_length, k, w, m, outbam=None):
    removed_reads = set()
    chrom_sizes = utils.get_chroms_list_from_bam(bam)
    samfile = pysam.AlignmentFile(bam, "rb") # type: ignore
    introns = get_introns_from_bam(samfile)
    index = index_samfile_by_reads(samfile)
    spurious_alignments, spurious_mates, spurious_introns, NH = find_spurious_alignments(introns, index, ref_fa,
                                                                        chrom_sizes, read_length, scoring, k, w, m)
    #write new bam
    if outbam is not None:
        removed_reads = write_filtered_bam(outbam, samfile, spurious_alignments, spurious_mates, NH, index)
    
    return spurious_introns, removed_reads



    

if __name__ == '__main__':
    import time
    scoring=[2,4,4,2,24,1,1]
    
    k = 7
    w = 7
    m = 22
    
    start = time.time()
    bam = "tests/data/ERR188044_chrX.bam"
    read_length = utils.get_read_length_from_bam(bam)
    ref_fa = "tests/data/chrX.fa"
    outbam = "tests/output/ERR188044_chrX_filtered.bam"

    spur_introns,removed_reads = filter_alignments_from_bam(ref_fa, bam, scoring, read_length, k, w,m)
    end = time.time()
    print(f"took {end-start} seconds")
    
    start = time.time()
    bam = "/ccb/salz8-2/shinder/projects/Geuvadis/BAM/ERR188025.bam"
    outbam = "/ccb/salz8-2/shinder/projects/Geuvadis/EASTR2/BAM/ERR188025_filtered.bam"
    out_introns = "/ccb/salz8-2/shinder/projects/Geuvadis/EASTR2/BED/ERR188025.bed"
    read_length = utils.get_read_length_from_bam(bam)
    outdir="/ccb/salz8-2/shinder/projects/Geuvadis/EASTR/BAM"
    ref_fa= "/ccb/salz8-2/chess-brain/ref/hg38mod_noPARs.fa"
    spurious_introns, removed_reads = filter_alignments_from_bam(ref_fa, bam, scoring, read_length, k, w, m, outbam=outbam)
    end = time.time()
    print(f"took {(end-start)/60} minutes")

    # #chroms = [x['SN'] for x in samfile.header['SQ']]
    
    # start = time.time()
    # bam = "/ccb/salz8-2/shinder/projects/repeat_spliced_alignments/chess3_tiebrush_bam/all_tissues.tb.bam"
    # introns = get_spurious_junctions_from_bam(bam, scoring, ref_fa, read_length,k,w)
    # end = time.time()
    # print(f"took {end-start} seconds")
    
    # introns[~introns['score'].isna()]