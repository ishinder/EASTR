import argparse
from EASTR import get_spurious_introns, utils, filter_bam
from io import StringIO
from posixpath import basename
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        prog="EASTR",
        description="Emend alignments of spuriously spliced transcript reads"
    )

    #required args
    parser.add_argument(
        "-R", "--reference",
        help="reference fasta genome used in alignment", required=True)
    
    parser.add_argument(
        "-bam",required=True,
        help="Input BAM file to emend alignments")

    #minimap2 args
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

    #output args
    parser.add_argument("--out_introns", default='stdout',metavar='FILE',
                        help="write introns to FILE; the default output is to terminal")

    parser.add_argument("--out_bam", metavar='FILE',
                        help="write filtered bam to FILE")

    #other args
    parser.add_argument(
        "-p",
        help="Number of parallel processes, default=1",
        default=1
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

def main():
    args = parse_args()
    # p = args.p
    scoring = minimap_scoring(args)
    ref_fa = args.reference
    bam = args.bam

    utils.index_bam(bam)
    utils.index_fasta(ref_fa)
    
    read_length = utils.get_read_length_from_bam(bam)

    # introns = get_spurious_introns.run_junctions(bam, scoring, ref_fa, 
    #                                              read_length, args.k, args.w)
    
    # filter_bam.filter_alignments(introns, ref_fa, bam, args.o)
    
    spurious_introns = get_spurious_introns.filter_alignments_from_bam(
            ref_fa, bam, scoring, read_length, args.k, args.w, outbam=args.out_bam)
    
    df = pd.Series(spurious_introns).reset_index()
    if args.out_introns == "stdout":
        out_introns = StringIO()
          
    else:
        out_introns=args.out_introns

    df.to_csv(out_introns,header=False,index=False, sep='\t')

    #TODO: do not make filtered bam
    #TODO: get spur introns only

if __name__ == '__main__':
    main()
    

