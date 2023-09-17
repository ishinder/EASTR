import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        prog="EASTR",
        description="EASTR: Emending alignments of spuriously spliced transcript reads."
                    "The script takes GTF, BED, or BAM files as input and processes "
                    "them using the provided reference genome and BowTie2 index. It "
                    "identifies spurious junctions and filters the input data "
                    "accordingly."
    )

    # Required Arguments
    group_reqd = parser.add_mutually_exclusive_group(required=True)
    group_reqd.add_argument("--gtf",
                            help="Input GTF file containing transcript annotations")
    group_reqd.add_argument("--bed",
                            help="Input BED file with intron coordinates")
    group_reqd.add_argument("--bam",
                            help="Input BAM file or a TXT file containing a list of BAM"
                            " files with read alignments")
    parser.add_argument("-r", "--reference",
                        required=True,
                        help="reference FASTA genome used in alignment")
    parser.add_argument("-i","--bowtie2_index",
                        required=True,
                        help="Path to Bowtie2 index for the reference genome")


    # bt2 Arguments:
    parser.add_argument(
        "--bt2_k",
        help="Minimum number of distinct alignments found by bowtie2 for a junction to"
        " be considered spurious. Default: 10",
        default=10,
        type=int
    )
    # EASTR Arguments
    parser.add_argument(
        "-o",
        help="Length of the overhang on either side of the splice junction. "
        "Default: 50",
        default=50,
        type=int)

    parser.add_argument(
       "--min_duplicate_exon_length",
        help="Minimum length of the duplicated exon. Default: 27",
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
            "if the flanking regions are similar (based on mappy scoring matrix). "
            "Default: 1",
        default=1,
        type=int)

    parser.add_argument(
        "--trusted_bed",
        help="Path to a BED file path with trusted junctions, which will "
        "not be removed by EASTR."
    )

    parser.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="Display additional information during BAM filtering, "
        "including the count of total spliced alignments and removed alignments")


    parser.add_argument( #TODO: directory instead of store_true
        "--removed_alignments_bam",
        default=False,
        action="store_true",
        help="Write removed alignments to a BAM file")

    # Minimap2 Arguments
    group_mm2 = parser.add_argument_group("Minimap2 parameters")

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
        help="Gap extension penalty. A gap of length k costs min(O1+k*E1, O2+k*E2). "
        "Default = [2, 1]",
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

    # Output Arguments
    group_out = parser.add_argument_group("Output")

    group_out.add_argument("--out_original_junctions",
                           default=None,
                           metavar="OUT",
                           help="Write original junctions to the OUT file or directory")

    group_out.add_argument("--out_removed_junctions",
                           default="stdout",
                           metavar="OUT",
                           help="Write removed junctions to OUT file or directory; "
                           "the default output is to terminal")

    group_out.add_argument("--out_filtered_bam",
                           metavar="OUT",
                           default=None,
                           help="Write filtered bams to OUT file or directory")

    group_out.add_argument("--filtered_bam_suffix",
                           metavar="STR",
                           default="_EASTR_filtered",
                           help="Suffix added to the name of the output BAM files. "
                           "Default=_EASTR_filtered")


    # Other Arguments
    parser.add_argument(
        "-p",
        help="Number of parallel processes, default=1",
        default=1,
        type=int
    )


    return parser.parse_args()
