import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        prog="EASTR",
        description="Emend alignments of spuriously spliced transcript reads"
    )

    parser.add_argument(
        "-R", "--reference",
        help="reference fasta genome used in alignment")
    
    parser.add_argument(
        "-B", "--bam",
        help="Input BAM file to emend alignments")

    return parser.parse_args()


def main():
    args = parse_args()
    REF_FA = args.reference
    BAM = args.bam
    print(args)
    print(REF_FA, BAM)


if __name__ == '__main__':
    main()

