import argparse
from collections import defaultdict
import csv
from ctypes import Union
import multiprocessing
import os
import sys
from EASTR import extract_junctions, utils, filter_bam, get_spurious_introns, alignment_utils
from io import StringIO
from posixpath import basename
import pandas as pd
from EASTR import filter_bam
from itertools import repeat
import pandas as pd
import re
from typing import List, Union

def out_junctions_filelist(bam_list:list, out_junctions, suffix="") -> Union[List[str], None]:

    if out_junctions is None:
        return

    if len(bam_list) == 1:
        if os.path.splitext(out_junctions)[1] != '.bed':
            out_junctions = out_junctions + '.bed'
        path = os.path.dirname(out_junctions)
        utils.make_dir(path)
        out_junctions_filelist = [out_junctions]

    else:
        if utils.check_directory_or_file(out_junctions) == 'file':
            print("ERROR: the path provided for the output file is a file, not a directory")
            sys.exit(1)
        utils.make_dir(out_junctions)
        out_junctions_filelist = []
        for bam in bam_list:
            out_junctions_filelist.append(out_junctions + "/" + os.path.splitext(os.path.basename(bam))[0] + suffix + ".bed")

    return out_junctions_filelist


def write_spurious_dict_to_all_stdout(spurious_dict, scoring):
    sorted_keys = spurious_dict.keys()
    named_keys = {}
    for i, key in enumerate(sorted_keys):
        name = "JUNC{}".format(i+1)
        named_keys[key] = name
    out = StringIO()    
    writer = csv.writer(out, delimiter='\t')
    for key, value in spurious_dict.items():
            chrom , start, end, strand = key
            name = named_keys[key]
            score = value['score']
            samples = value['samples']
            score2 = alignment_utils.calc_alignment_score(value['hit'],scoring)
            name2 = ';'.join([f"{score2},{sample_id}" for score2, sample_id in samples])
            writer.writerow([chrom, start, end, name, score, strand, score2, name2])
    print(out.getvalue())


def write_spurious_dict_to_sample_bed(spurious_dict, bam_list, out_removed_junctions_filelist, scoring):
    sorted_keys = spurious_dict.keys()
    named_keys = {}
    for i, key in enumerate(sorted_keys):
        name = "JUNC{}".format(i+1)
        named_keys[key] = name

    sample_names = [os.path.splitext(os.path.basename(bam_path))[0] for bam_path in bam_list]

    out_io_dict = {}
    for i,sample in enumerate(sample_names):
        out_io_dict[sample] = [StringIO(), out_removed_junctions_filelist[i]]

    for key, value in spurious_dict.items():
        name = named_keys[key]
        samples = value['samples']
        score2 = alignment_utils.calc_alignment_score(value['hit'],scoring)
        for i, sample in enumerate(samples):
            score, sample_id = sample
            out_io_dict[sample_id][0].write(f"{key[0]}\t{key[1]}\t{key[2]}\t{name}\t{score}\t{key[3]}\t{score2}")

    for (out,filepath) in out_io_dict.values():
        with open(filepath, 'w') as out_removed_junctions:
            out_removed_junctions.write(out.getvalue())
        


def parse_args(arglist):
    parser = argparse.ArgumentParser(
        prog="EASTR",
        description="Emending alignments of spuriously spliced transcript reads"
    )

    #required args
    group_reqd = parser.add_mutually_exclusive_group(required=True)
    group_reqd.add_argument('--gtf', help='Input GTF file')
    group_reqd.add_argument('--bed', help='Input BED file with intron coodrdiantes')
    group_reqd.add_argument('--bam', help='Input BAM file or a TXT file containing a list of BAM files')
    parser.add_argument("-r", "--reference", required=True, help="reference fasta genome used in alignment")
    parser.add_argument('-i','--bowtie2_index', required=True, help='path to bowtie2 index')

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
        "--min_junc_score", 
        help=" minimum number of supporting spliced reads required per junction. \
        Any self-aligning junction with less than the minimum number of supporting reads per all samples is filtered",
        default=1,
        type=int)

    parser.add_argument(
        "--trusted_bed", 
        help="BED file path with trusted junctions (will not be removed by EASTR)")
    


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

    group_out.add_argument("--out_original_junctions", default=None, metavar='OUT',
                        help="write original junctions to OUT file or directory")
    
    group_out.add_argument("--out_removed_junctions", default='stdout', metavar='FILE',
                        help="write removed junctions to OUT file or directory; the default output is to terminal")

    group_out.add_argument("--out_filtered_bam", metavar='OUT',
                        help="write filtered bams to OUT file or directory")

    group_out.add_argument("--filtered_bam_suffix", metavar='STR',default="_EASTR_filtered",
                        help="suffix added to the name of the output BAM files, default='_EASTR_filtered'")


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
    
    #required input args
    bam_list = args.bam
    gtf_path = args.gtf
    bed_path = args.bed
    ref_fa = args.reference
    bt2_index = args.bowtie2_index
    bt2_k = args.bt2_k

    #EASTR variables
    overhang = args.o
    min_junc_score = args.min_junc_score
    anchor = args.a
    trusted_bed = args.trusted_bed
    
    #mm2 variables
    scoring = minimap_scoring(args)
    k = args.k
    w = args.w
    m = args.m

    #output args
    suffix = args.filtered_bam_suffix
    out_original_junctions = args.out_original_junctions
    out_removed_junctions = args.out_removed_junctions
    out_filtered_bam = args.out_filtered_bam


    #other args
    p = args.p

    #index reference fasta if it's not indexed
    utils.index_fasta(ref_fa)

    #check if the input is a bam file or a list of bam files
    is_bam = False
    if bam_list:
        is_bam = True

        #if a single bam file is provided
        extension = os.path.splitext(os.path.basename(bam_list))[1]
        if extension in ['.bam','.cram','.sam']:
            bam_list = [bam_list]
        
        #if a file containing a list of bam files is provided
        else:
            with open(bam_list) as file:
                bam_list = [line.rstrip() for line in file]
                

    out_original_junctions_filelist = out_junctions_filelist(bam_list, out_original_junctions, suffix=suffix)
    out_removed_junctions_filelist = out_junctions_filelist(bam_list, out_removed_junctions, suffix=suffix)
    
    spurious_dict = get_spurious_introns.get_spurious_junctions(scoring, k, w, m, overhang, bt2_index, bt2_k,
                                    ref_fa, p, anchor, min_junc_score, bam_list, gtf_path,
                                    bed_path, trusted_bed, out_original_junctions_filelist )




    if out_removed_junctions_filelist is None:
        out_removed_junctions = StringIO()
        #TODO
          
   

    df = get_spurious_introns.spurious_junctions_to_bed(spurious_dict, out_removed_junctions, is_bam = is_bam, scoring=scoring)

    if bam_list:

        if args.out_dir:
            out_dir = args.out_dir
            utils.make_dir(out_dir)
            outbams = [os.path.join(out_dir, re.split(r'.[bcr]+am',basename(i))[0] + f'{suffix}.bam') for i in bam_list]

            pool = multiprocessing.Pool(p)
            removed_reads = pool.starmap(filter_bam.write_filtered_bam, zip(bam_list,outbams,repeat(set(spurious_dict))))
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
    

