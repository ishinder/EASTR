import argparse
from collections import defaultdict
import csv
import multiprocessing
import os
from EASTR import extract_junctions, utils, filter_bam, get_spurious_introns
from io import StringIO
from posixpath import basename
import pandas as pd
from EASTR import filter_bam
from itertools import repeat
import pandas as pd
import re

def parse_args(arglist):
    parser = argparse.ArgumentParser(
        prog="EASTR",
        description="Emending alignments of spuriously spliced transcript reads"
    )

    #required args
    group_reqd = parser.add_mutually_exclusive_group(required=True)
    group_reqd.add_argument('--gtf', help='Input GTF file to identify potentially spurious introns')
    group_reqd.add_argument('--bed', help='Input BED file with intron coodrdiantes to identify potentially spurious introns')
    group_reqd.add_argument('--bam', help='Input BAM file to identify potentially spurious spliced alignments')

    parser.add_argument("-r", "--reference", required=True, help="reference fasta genome used in alignment")
    parser.add_argument( '-i','--bowtie2_index', required=True, help='path to bowtie2 index')

    #gtf args:
    parser.add_argument("--trust_gtf", action='store_true')

    #bt2 args:
    parser.add_argument(
        "--bt2_k",
        help="minimum number of distinct alignments found by bowtie2 such that a given junction may be \
        considered spurious",
        default=10,
        type=int
    )
    #EASTR args
    parser.add_argument(
        "-o",
        help="overhang on either side of the splice junction, default = 50",
        default=50,
        type=int)

    parser.add_argument(
        "-a",
        help="minimum required anchor length in each of the two exons, default = 7",
        default=7,
        type=int
    )

    parser.add_argument(
        "--min_junc_score_global", 
        help=" minimum number of supporting spliced reads required per junction. \
        Any self-aligning junction with less than the minimum number of supporting reads per all samples is filtered",
        default=1,
        type=int)

    parser.add_argument(
        "--min_junc_score_local", 
        help=" minimum number of supporting reads required per junction. \
        Any self-aligning junction with less than the minimum number of supporting reads per sample is filtered",
        default=1,
        type=int)

    parser.add_argument(
        "--trusted_gtf", 
        help="GTF file path with trusted junctions (will not be filtered)")
    


    #minimap2 args
    group_mm2 = parser.add_argument_group('Minimap2 parameters')

    group_mm2.add_argument(
        "-A",
        help="Matching score, default = 3",
        default=3,
        type=int)

    group_mm2.add_argument(
        "-B",
        help="Mismatching penalty, default = 4",
        default=4,
        type=int)

    group_mm2.add_argument(
        "-O",
        nargs=2,
        type=int,
        help="Gap open penalty, default = [12, 32]",
        default=[12,32])
    
    group_mm2.add_argument(
        "-E",
        nargs=2,
        type=int,
        help="Gap extension penalty, default = [2, 1]. A gap of length k costs min(O1+k*E1, O2+k*E2).",
        default=[2,1])

    group_mm2.add_argument(
        "-k",
        help="kmer length for alignment, default=3",
        default=3,
        type=int
    )

    group_mm2.add_argument(
        "--scoreN",
        help="Score of a mismatch involving ambiguous bases, default=1",
        default=1,
        type=int
    )

    group_mm2.add_argument(
        "-w",
        help="minimizer window size, default=2",
        default=2, 
        type=int
    )

    group_mm2.add_argument(
        "-m",
        help="Discard chains with chaining score, default=25.",
        default=25, 
        type=int
    )

    #output args
    group_out = parser.add_argument_group('Output')

    group_out.add_argument("--out_introns", default='stdout',metavar='FILE',
                        help="write spurious introns to FILE; the default output is to terminal")

    group_out.add_argument("--out_dir", metavar='DIRECTORY',
                        help="write filtered bams and removed reads to DIRECTORY")

    group_out.add_argument("--suffix", metavar='STR',default="_filtered",
                        help="suffix added to output BAM files, default='.suffix'")


    #other args
    parser.add_argument(
        "-p",
        help="Number of parallel processes, default=1",
        default=1,
        type=int
    )



    return parser.parse_args()

def minimap_scoring(args):
    gap_open_penalty = args.O
    gap_ext_penalty = args.E
    mismatch_penalty = args.B
    match_score = args.A
    ambiguous_score = args.scoreN

    scoring=[
        match_score,
        mismatch_penalty,
        gap_open_penalty[0], 
        gap_ext_penalty[0],
        gap_open_penalty[1], 
        gap_ext_penalty[1],
        ambiguous_score]
    
    return scoring

#def write_removed_reads(outfile, removed_reads):


def main(arglist=None):
    args = parse_args(arglist)
    
    #set values
    scoring = minimap_scoring(args)
    k = args.k
    w = args.w
    m = args.m
    overhang = args.o
    min_junc_score_global = args.min_junc_score_global
    min_junc_score_local = args.min_junc_score_local
    anchor = args.a
    p = args.p
    ref_fa = args.reference
    bt2_index = args.bowtie2_index
    bt2_k = args.bt2_k
    suffix = args.suffix

    is_bam = False
    trusted_introns = args.trusted_gtf
    bam_list = args.bam
    gtf_path = args.gtf
    bed_path = args.bed

    utils.index_fasta(ref_fa)

    if bam_list is not None:
        is_bam = True

        #if a single bam file is provided
        if bam_list.split('.')[-1] in ['bam','cram']:
            bam_list = [bam_list]
        
        #if a file containing a list of bam files is provided
        else:
            with open(bam_list) as file:
                bam_list = [line.rstrip() for line in file]
                
        pool = multiprocessing.Pool(p)
        d = pool.map(utils.index_bam,bam_list)
    
    
    spurious = get_spurious_introns.get_spurious_junctions(scoring, k, w, m, overhang, bt2_index, bt2_k,
                                    ref_fa, p, anchor, min_junc_score_global, min_junc_score_local, bam_list=bam_list, 
                                    gtf_path=gtf_path, trusted_introns=trusted_introns)

    
    print(args.out_introns)
    if args.out_introns == "stdout":
            out_introns = StringIO()
          
    else:
        out_introns = args.out_introns
        path = os.path.dirname(out_introns)
        if path != '':
            utils.make_dir(os.path.dirname(out_introns))

    df = get_spurious_introns.spurious_junctions_to_bed(spurious, out_introns, is_bam = is_bam, scoring=scoring)

    if bam_list:

        if args.out_dir:
            out_dir = args.out_dir
            utils.make_dir(out_dir)
            outbams = [os.path.join(out_dir, re.split(r'.[bcr]+am',basename(i))[0] + f'{suffix}.bam') for i in bam_list]

            pool = multiprocessing.Pool(p)
            removed_reads = pool.starmap(filter_bam.write_filtered_bam, zip(bam_list,outbams,repeat(set(spurious))))
            pool.close()

            for i, rr in enumerate(removed_reads):
                out_removed_reads = os.path.join(out_dir, re.split(r'.[bcr]+am',basename(bam_list[i]))[0] + '.removed_reads.lst')
                pd.DataFrame(rr).to_csv(out_removed_reads, index=False,header=False,sep='\t')

            for i, bam in enumerate(bam_list):
                sample_id = re.split(r'.[bcr]+am',basename(bam))[0]
                
                #write removed_reads
                out_removed_reads = os.path.join(out_dir, sample_id + '.removed_reads.lst')
                pd.DataFrame(removed_reads[i]).to_csv(out_removed_reads, index=False, header=False, sep='\t')
                
                #write bed file of spurious junctions
                df_sample = df[df['samples'].str.contains(sample_id)].iloc[:,:-1]
                out_spur_juncs = os.path.join(out_dir, sample_id + '.spurious_junctions.bed')
                df_sample.to_csv(out_spur_juncs, header=False, index=False, sep='\t')

if __name__ == '__main__':
    main()
    

