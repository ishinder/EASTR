import argparse
import os

from EASTR import get_spurious_introns
from EASTR import output_utils
from EASTR import utils

def parse_args(arglist):
    parser = argparse.ArgumentParser(
        prog="EASTR",
        description="EASTR: Emending alignments of spuriously spliced transcript reads. "
                    "The script takes GTF, BED, or BAM files as input and processes them using "
                    "the provided reference genome and BowTie2 index. It identifies spurious junctions "
                    "and filters the input data accordingly."
    )

    #required args
    group_reqd = parser.add_mutually_exclusive_group(required=True)
    group_reqd.add_argument('--gtf', help='Input GTF file containing transcript annotations')
    group_reqd.add_argument('--bed', help='Input BED file with intron coordinates')
    group_reqd.add_argument('--bam', help='Input BAM file or a TXT file containing a list of BAM files with read alignments')
    parser.add_argument("-r", "--reference", required=True, help="reference FASTA genome used in alignment")
    parser.add_argument('-i','--bowtie2_index', required=True, help='Path to Bowtie2 index for the reference genome')


    #bt2 args:
    parser.add_argument(
        "--bt2_k",
        help="Minimum number of distinct alignments found by bowtie2 for a junction to be \
        considered spurious. Default: 10",
        default=10,
        type=int
    )
    #EASTR args
    parser.add_argument(
        "-o",
        help="Length of the overhang on either side of the splice junction. Default = 50",
        default=50,
        type=int)

    parser.add_argument(
       "--min_duplicate_exon_length",
         help="Minimum length of the duplicated exon. Default = 27",
        default=27,
        type=int
    )

    parser.add_argument(
        "-a",
        help="Minimum required anchor length in each of the two exons, default = 7",
        default=7,
        type=int
    )

    parser.add_argument(
        "--min_junc_score",
        help=" Minimum number of supporting spliced reads required per junction. "
            "Junctions with fewer supporting reads in all samples are filtered out "
            "if the flanking regions are similar (based on mappy scoring matrix). Default: 1",
        default=1,
        type=int)

    parser.add_argument(
        "--trusted_bed",
        help="Path to a BED file path with trusted junctions, which will not be removed by EASTR."
    )

    parser.add_argument(
        "--verbose", default=False, action="store_true",
        help="Display additional information during BAM filtering, "
        "including the count of total spliced alignments and removed alignments")


    parser.add_argument( #TODO: directory instead of store_true
        "--removed_alignments_bam", default=False, action="store_true",
        help="Write removed alignments to a BAM file")

    #minimap2 args
    group_mm2 = parser.add_argument_group('Minimap2 parameters')

    group_mm2.add_argument(
        "-A",
        help="Matching score. Default = 3",
        default=3,
        type=int)

    group_mm2.add_argument(
        "-B",
        help="Mismatching penalty. Default = 4",
        default=4,
        type=int)

    group_mm2.add_argument(
        "-O",
        nargs=2,
        type=int,
        help="Gap open penalty. Default = [12, 32]",
        default=[12,32])

    group_mm2.add_argument(
        "-E",
        nargs=2,
        type=int,
        help="Gap extension penalty. A gap of length k costs min(O1+k*E1, O2+k*E2). Default = [2, 1]",
        default=[2,1])

    group_mm2.add_argument(
        "-k",
        help="K-mer length for alignment. Default=3",
        default=3,
        type=int
    )

    group_mm2.add_argument(
        "--scoreN",
        help="Score of a mismatch involving ambiguous bases. Default=1",
        default=1,
        type=int
    )

    group_mm2.add_argument(
        "-w",
        help="Minimizer window size. Default=2",
        default=2,
        type=int
    )

    group_mm2.add_argument(
        "-m",
        help="Discard chains with chaining score. Default=25.",
        default=25,
        type=int
    )

    #output args
    group_out = parser.add_argument_group('Output')

    group_out.add_argument("--out_original_junctions", default=None, metavar='OUT',
                        help="Write original junctions to the OUT file or directory")

    group_out.add_argument("--out_removed_junctions", default='stdout', metavar='OUT',
                        help="Write removed junctions to OUT file or directory; the default output is to terminal")

    group_out.add_argument("--out_filtered_bam", metavar='OUT',default=None,
                        help="Write filtered bams to OUT file or directory")

    group_out.add_argument("--filtered_bam_suffix", metavar='STR',default="_EASTR_filtered",
                        help="Suffix added to the name of the output BAM files. Default='_EASTR_filtered'")


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
    bed_list = args.bed
    ref_fa = args.reference
    bt2_index = args.bowtie2_index
    bt2_k = args.bt2_k

    #EASTR variables
    overhang = args.o
    min_duplicate_exon_length = args.min_duplicate_exon_length
    min_junc_score = args.min_junc_score
    anchor = args.a
    trusted_bed = args.trusted_bed
    verbose = args.verbose
    removed_alignments_bam = args.removed_alignments_bam

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

        else:
            with open(bam_list) as file:
                bam_list = [line.rstrip() for line in file]
                for bam in bam_list:
                    if not os.path.isfile(bam):
                        raise ValueError('input must be a bam file or a file containing a list of bam files')


    elif bed_list:
        extension = os.path.splitext(os.path.basename(bed_list))[1]
        if extension in ['.bed']:
            bed_list = [bed_list]

        else:
            with open(bed_list) as file:
                bed_list = [line.rstrip() for line in file]
                for bed in bed_list:
                    if not os.path.isfile(bed):
                        raise ValueError('input must be a bed file or a file containing a list of bed files')


    original_junctions_filelist = output_utils.out_junctions_filelist(bam_list, gtf_path, bed_list, out_original_junctions, suffix="_original_junctions")
    removed_junctions_filelist = output_utils.out_junctions_filelist(bam_list, gtf_path, bed_list, out_removed_junctions, suffix="_removed_junctions")
    filtered_bam_filelist = output_utils.out_filtered_bam_filelist(bam_list, out_filtered_bam, suffix=suffix)

    spurious_dict = get_spurious_introns.get_spurious_junctions(scoring, k, w, m, overhang, min_duplicate_exon_length, bt2_index, bt2_k,
                                    ref_fa, p, anchor, min_junc_score, bam_list, gtf_path,
                                    bed_list, trusted_bed, original_junctions_filelist, verbose )


    if is_bam:
        if filtered_bam_filelist:
            sample_to_bed = output_utils.spurious_dict_bam_by_sample_to_bed(
                    spurious_dict, bam_list, removed_junctions_filelist, scoring=scoring)
            output_utils.filter_multi_bam_with_vacuum(bam_list, sample_to_bed, filtered_bam_filelist, p, verbose, removed_alignments_bam)
            if removed_junctions_filelist is None:
                for _, sample in sample_to_bed.items():
                    os.remove(sample)
        elif removed_junctions_filelist:
            sample_to_bed = output_utils.spurious_dict_bam_by_sample_to_bed(spurious_dict, bam_list, removed_junctions_filelist, scoring=scoring)

        else:
            output_utils.spurious_dict_all_to_bed(spurious_dict, scoring, None, gtf_path, bed_list, bam_list)

    if gtf_path:
        output_utils.spurious_dict_all_to_bed(spurious_dict, scoring, removed_junctions_filelist, gtf_path, bed_list, bam_list)

    elif bed_list:
        if removed_junctions_filelist:
            output_utils.spurious_dict_bed_by_sample_to_bed(spurious_dict, bed_list, removed_junctions_filelist, scoring)
        else:
            output_utils.spurious_dict_all_to_bed(spurious_dict, scoring, None, gtf_path, bed_list, bam_list)

if __name__ == '__main__':
    main()
