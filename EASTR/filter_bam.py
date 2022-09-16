import pysam
import pandas as pd
import mappy as mp
import collections
from posixpath import basename, splitext, dirname
import os


def make_sam_header(ref_fa,VN='1.0'):
    fasta = pysam.FastaFile(ref_fa) # type: ignore
    chroms = fasta.references
    header = {}
    header['HD'] = {'VN':VN}
    header['SQ'] = []

    for SN in chroms:
        LN = fasta.get_reference_length(SN)
        header['SQ'].append({'LN':LN, 'SN':SN})

    return header

def write_list_to_file(input, outfile):
    with open(outfile, mode='wt', encoding='utf-8') as f:
        f.write('\n'.join(input))
        f.write('\n')
    f.close()
    
#TODO modify NH tag for kept alignments

def filter_alignments(introns, ref_fa, inbam, outdir):   
    name = basename(inbam)
    name = ''.join(name.split('.')[:-1])
    outbam_filtered = outdir + "/" + name + "_filtered.bam"
    outlist_reads = outdir + "/" + name + "_removed_reads.lst"
    
    #get introns that may be outliers
    min_score = min(introns['max_o_score'])
    filt_introns = introns[introns['max_o_score']!=min_score]
    
    header = make_sam_header(ref_fa)
    spurAlignments = set() 
    keepAlignments = collections.defaultdict(list)
    samfile = pysam.AlignmentFile(inbam, "rb") # type: ignore

    for alignment in samfile.fetch(until_eof=True):
        if len(alignment.blocks) == 1:
            keepAlignments[alignment.qname].append(alignment)
            continue
        
        keep = True
        currentloc = alignment.pos
        for i,cigarop in enumerate(alignment.cigar):
            if (cigarop[0]==4): #substitution or insertion in query
                continue
            if (cigarop[0]==1):
                continue
            if(cigarop[0]==3):
                key = (samfile.getrname(alignment.tid),currentloc,currentloc+cigarop[1])
                if key in filt_introns.index:
                    keep = False
                    if not alignment.is_secondary:
                        spurAlignments.add(alignment.qname)
                        
            currentloc=currentloc+cigarop[1]
            
        if keep:
            keepAlignments[alignment.qname].append(alignment)     
        
    writeKeys = keepAlignments.keys() - spurAlignments
    writeAlignments = {key: keepAlignments[key] for key in writeKeys}

    #make new bam
    outf = pysam.AlignmentFile(outbam_filtered, "wb", header=header) # type: ignore
    for alignments in writeAlignments.values():
        for alignment in alignments:
            w=outf.write(alignment)
    outf.close()
    x=pysam.sort("-o", outbam_filtered, outbam_filtered)
    
    write_list_to_file(spurAlignments, outlist_reads)
    

if __name__ == '__main__':
    
    import time
    from EASTR import get_spurious_introns, utils
    
    start = time.time()
    outdir = 'tests/output/'
    bam = "tests/data/ERR188044_chrX.bam"
    ref_fa = "tests/data/chrX.fa"
    scoring=[2,4,4,2,24,1,1]
    read_length = utils.get_read_length_from_bam(bam)
    k = 7
    w = 7
    introns = get_spurious_introns.run_junctions(bam, scoring, ref_fa, read_length,k,w)
    filter_alignments(introns, ref_fa, bam, outdir)
    end = time.time()
    print(f"took {end-start} seconds") #107 seconds

