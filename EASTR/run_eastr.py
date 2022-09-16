import argparse
from EASTR import get_spurious_introns, utils, filter_bam
from io import StringIO
from posixpath import basename

def parse_args():
    parser = argparse.ArgumentParser(
        prog="EASTR",
        description="Emend alignments of spuriously spliced transcript reads"
    )

    parser.add_argument(
        "-R", "--reference",
        help="reference fasta genome used in alignment")
    
    parser.add_argument(
        "-bam",
        help="Input BAM file to emend alignments")

    parser.add_argument(
        "-A",
        help="Matching score, default = 2",
        default=2)

    parser.add_argument(
        "-B",
        help="Mismatching penalty, default = 4",
        default=2)

    parser.add_argument(
        "-O",
        nargs=2,
        type=int,
        help="Gap open penalty, default = [4, 24]",
        default=[4,24])
    
    parser.add_argument(
        "-E",
        nargs=2,
        type=int,
        help="Gap extension penalty, default = [2, 1]. A gap of length k costs min(O1+k*E1, O2+k*E2).",
        default=[2,1])

    parser.add_argument(
        "-k",
        help="kmer length for alignment, default=7",
        default=7
    )

    parser.add_argument(
        "--scoreN",
        help="Score of a mismatch involving ambiguous bases, default=1",
        default=1
    )

    parser.add_argument(
        "-w",
        help="minimizer window size, default=7",
        default=7
    )

    parser.add_argument(
        "-p",
        help="Number of parallel processes, default=1",
        default=1
    )
    
    parser.add_argument("-o", default='stdout',metavar='FILE',
                        help="write output to FILE; the default output is to terminal")


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

def main():
    args = parse_args()
    # p = args.p
    scoring = minimap_scoring(args)
    ref_fa = args.reference
    bam = args.bam

    utils.index_bam(bam)
    utils.index_fasta(ref_fa)
    
    read_length = utils.get_read_length_from_bam(bam)

    introns = get_spurious_introns.run_junctions(bam, scoring, ref_fa, 
                                                 read_length, args.k, args.w)
    filter_bam.filter_alignments(introns, ref_fa, bam, args.o)
    
    if args.o=="stdout": #TODO - this should be required. 
        output = StringIO()
        introns.to_csv(output)
        
    else:
        name = basename(bam)
        name = ''.join(name.split('.')[:-1])
        introns.to_csv(args.o + f"/{name}_junctions.tsv",sep='\t')


if __name__ == '__main__':
    main()
    

