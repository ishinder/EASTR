import os
from collections import namedtuple

from EASTR import get_spurious_introns, output_utils, utils
from EASTR.helpers.args import parse_args

input_args = ["args"
        "bam_list"
        "bed_list"
        "filtered_bam_filelist"
        "removed_junctions_filelist"
        "scoring"
        "spurious_dict"]
Eastr_args = namedtuple("Eastr_args", input_args)

def handle_bam(e: Eastr_args):
    if e.filtered_bam_filelist:
        sample_to_bed = output_utils.spurious_dict_bam_by_sample_to_bed(
                    bam_list=e.bam_list,
                    out_removed_junctions_filelist=e.removed_junctions_filelist,
                    scoring=e.scoring,
                    spurious_dict=e.spurious_dict,)
        output_utils.filter_multi_bam_with_vacuum(bam_list=e.bam_list,
                                                  sample_to_bed=sample_to_bed,
                                                  out_bam_list=e.filtered_bam_filelist,
                                                  num_parallel=e.args.p,
                                                  verbose=e.args.verbose,
                                                  removed_alignments_bam=e.args.removed_alignments_bam)
        if e.removed_junctions_filelist is None:
            for sample in sample_to_bed:
                os.remove(sample_to_bed[sample])
    elif e.removed_junctions_filelist:
        sample_to_bed = output_utils.spurious_dict_bam_by_sample_to_bed(e.spurious_dict,
                                                                        e.bam_list,
                                                                        e.removed_junctions_filelist,
                                                                        scoring=e.scoring)
    else:
        output_utils.spurious_dict_all_to_bed(bam_list=e.bam_list,
                                              bed_list=e.bed_list,
                                              gtf_path=e.args.gtf,
                                              fileout=None,
                                              scoring=e.scoring,
                                              spurious_dict=e.spurious_dict,)


def handle_bed_list(e: Eastr_args):
    if e.removed_junctions_filelist:
        output_utils.spurious_dict_bed_by_sample_to_bed(bed_list=e.bed_list,
                                                        out_removed_junctions_filelist=e.removed_junctions_filelist,
                                                        scoring=e.scoring,
                                                        spurious_dict=e.spurious_dict)
    else:
        output_utils.spurious_dict_all_to_bed(bam_list=e.bam_list,
                                              bed_list=e.bed_list,
                                              fileout=None,
                                              gtf_path=e.args.gtf_path,
                                              scoring=e.scoring,
                                              spurious_dict=e.spurious_dict)


def handle_gtf(e: Eastr_args):
    output_utils.spurious_dict_all_to_bed(bam_list=e.bam_list,
                                          bed_list=e.bed_list,
                                          fileout=e.removed_junctions_filelist,
                                          gtf_path=e.args.gtf_path,
                                          scoring=e.scoring,
                                          spurious_dict=e.spurious_dict)


def main():
    args = parse_args()

    # Required input args
    bam_list = args.bam
    gtf_path = args.gtf
    bed_list = args.bed


    # mm2 variables
    scoring = utils.get_minimap_scoring_list(args)

    # output args
    out_original_junctions = args.out_original_junctions

    # index reference fasta if it's not indexed
    utils.index_fasta(args.reference)

    # check if the input is a bam file or a list of bam files
    if bam_list:
        bam_list = utils.get_bam_list(bam_list)
    elif bed_list:
        bed_list = utils.get_bed_list(bed_list)

    original_junctions_filelist = output_utils.out_junctions_filelist(bam_list,
                                                                      gtf_path,
                                                                      bed_list,
                                                                      out_original_junctions,
                                                                      suffix="_original_junctions")
    removed_junctions_filelist = output_utils.out_junctions_filelist(bam_list,
                                                                     gtf_path,
                                                                     bed_list,
                                                                     args.out_removed_junctions,
                                                                     suffix="_removed_junctions")
    filtered_bam_filelist = output_utils.out_filtered_bam_filelist(bam_list,
                                                                   args.out_filtered_bam,
                                                                   suffix=args.filtered_bam_suffix)

    spurious_dict = get_spurious_introns.get_spurious_junctions(args=args,
                                                                scoring=scoring,
                                                                bam_list=bam_list,
                                                                bed_path=bed_list,
                                                                out_original_junctions=original_junctions_filelist
                                                                )

    eastr_input_args = Eastr_args(args=args,
                                  bam_list=bam_list,
                                  bed_list=bed_list,
                                  filtered_bam_filelist=filtered_bam_filelist,
                                  removed_junctions_filelist=removed_junctions_filelist,
                                  scoring=scoring,
                                  spurious_dict=spurious_dict)
    if bam_list:
        handle_bam(eastr_input_args)

    elif gtf_path:
        handle_gtf(eastr_input_args)

    elif bed_list:
        handle_bed_list(eastr_input_args)


if __name__ == "__main__":
    main()

