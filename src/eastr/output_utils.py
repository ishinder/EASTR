import io
import csv
import multiprocessing
import os
import pathlib
import shlex
import subprocess
import sys
import tempfile
from typing import List, Optional, Union

from eastr import alignment_utils
from eastr import utils

this_directory = pathlib.Path(__file__).resolve().parent
# This should exist with source after compilation.
VACUUM_CMD = os.path.join(this_directory, 'vacuum')

def out_junctions_filelist(
        bam_list: Optional[list[str]] = None,
        bed_list: Optional[list[str]] = None,
        gtf_path: Optional[str] = None,
        out_junctions: Optional[str] = None,
        suffix: Optional[str] = None,
) -> Union[List[str], None, str]:
    if out_junctions is None:
        return None

    if gtf_path:
        if os.path.isdir(out_junctions):
            out_junctions = os.path.join(
                out_junctions,
                utils.sanitize_and_update_extension(gtf_path, f"{suffix}.bed"),
            )
        return out_junctions

    if bed_list:
        if suffix in {"_original_junctions", "", None}:
            return None

        if len(bed_list) == 1:
            if os.path.isdir(out_junctions):
                out_junctions = os.path.join(
                    out_junctions,
                    utils.sanitize_and_update_extension(gtf_path, f"{suffix}.bed"),
                )
            path = os.path.dirname(out_junctions)
            utils.make_dir(path)
            return [out_junctions]

        if not os.path.isdir(out_junctions):
            print("ERROR: the path provided for the output bed files is a file path, not a directory")
            sys.exit(1)

        utils.make_dir(out_junctions)
        result = []
        for bed in bed_list:
            result.append(
                os.path.join(
                    out_junctions,
                    utils.sanitize_and_update_extension(bed, f"{suffix}.bed"),
                )
            )
        return result

    if len(bam_list) == 1:
        if os.path.isdir(out_junctions):
            out_junctions = os.path.join(
                out_junctions,
                utils.sanitize_and_update_extension(bam_list[0], f"{suffix}.bed"),
            )
        path = os.path.dirname(out_junctions)
        utils.make_dir(path)
        return [out_junctions]

    if not os.path.isdir(out_junctions):
        print("ERROR: the path provided for the output bed files is a file path, not a directory")
        sys.exit(1)

    utils.make_dir(out_junctions)
    result = []
    for bam in bam_list:
        result.append(
            os.path.join(
                out_junctions,
                utils.sanitize_and_update_extension(bam, f"{suffix}.bed"),
            )
        )
    return result

def out_filtered_bam_filelist(
        bam_list: Optional[list[str]] = None,
        out_filtered_bam: Optional[str] = None,
        suffix: Optional[str] = None,
) -> Union[List[str], None]:
    result = None
    if bam_list is None or out_filtered_bam is None:
        return result

    if suffix is None:
        suffix = "_EASTR_filtered"

    if len(bam_list) == 1:
        if os.path.isdir(out_filtered_bam):
            out_filtered_bam = os.path.join(
                out_filtered_bam,
                utils.sanitize_and_update_extension(bam_list[0], f"{suffix}.bam"),
            )
        path = os.path.dirname(out_filtered_bam)
        utils.make_dir(path)
        result = [out_filtered_bam]
    else:
        if not os.path.isdir(out_filtered_bam):
            print("ERROR: the path provided for the output file is a file, not a directory")
            sys.exit(1)
        utils.make_dir(out_filtered_bam)
        result = []
        for bam in bam_list:
            result.append(
                os.path.join(
                out_filtered_bam,
                utils.sanitize_and_update_extension(bam, f"{suffix}.bam"),
                )
            )
    return result

def writer_spurious_dict_bam_to_bed(spurious_dict, named_keys, scoring, writer):
    for key, value in spurious_dict.items():
        chrom, start, end, strand = key
        name = named_keys[key]
        score = value['score']
        samples = value['samples']
        score2 = alignment_utils.calc_alignment_score(value['hit'],scoring)
        name2 = ';'.join([f"{score2},{sample_id}" for score2, sample_id in samples])
        writer.writerow([chrom, start, end, name, score, strand, score2, name2])


def writer_spurious_dict_gtf_to_bed(spurious_dict, named_keys, scoring, writer):
    for key, value in spurious_dict.items():
        chrom, start, end, strand = key
        gene_id = value['transcripts'][0]
        transcripts = value['transcripts'][1]
        name = named_keys[key]
        score = '.'
        score2 = alignment_utils.calc_alignment_score(value['hit'],scoring)
        name2 = ';'.join([f"{gene_id}", *transcripts])
        writer.writerow([chrom, start, end, name, score, strand, score2, name2])

def writer_spurious_dict_bed_to_bed(spurious_dict, named_keys, scoring, writer):
    for key, value in spurious_dict.items():
        chrom, start, end, strand = key
        name = named_keys[key]
        samples = value['samples']
        score = value['score']
        name2 = ';'.join([f"{name2},{score2}" for _, name2, score2 in samples])
        writer.writerow([chrom, start, end, name, score, strand, name2])


def spurious_dict_all_to_bed(spurious_dict,scoring,fileout,gtf_path, bed_list, bam_list):
    sorted_keys = spurious_dict.keys()
    named_keys = {}
    for i, key in enumerate(sorted_keys):
        name = "JUNC{}".format(i+1)
        named_keys[key] = name
    out = io.StringIO()
    writer = csv.writer(out, delimiter='\t')

    if gtf_path is not None:
        writer_spurious_dict_gtf_to_bed(spurious_dict, named_keys, scoring, writer)
    elif bed_list is not None:
        writer_spurious_dict_bed_to_bed(spurious_dict, named_keys, scoring, writer)
    else:
        writer_spurious_dict_bam_to_bed(spurious_dict, named_keys, scoring, writer)

    if fileout is None:
        print(out.getvalue())
    else:
        with open(fileout, 'w') as out_file:
            out_file.write(out.getvalue())


def create_sample_to_bed_dict(sample_names, out_removed_junctions_filelist):
    sample_to_bed = {}

    for sample in sample_names:
        shortest_match = None
        for sample_id in out_removed_junctions_filelist:
            if sample in sample_id:
                if shortest_match is None or len(sample_id) < len(shortest_match):
                    shortest_match = sample_id
        if shortest_match is not None:  # Only update if a match was found
            file_path = out_removed_junctions_filelist[out_removed_junctions_filelist.index(shortest_match)]
            file_obj = open(file_path, mode='w+b')
            sample_to_bed[sample] = file_obj

    return sample_to_bed



def spurious_dict_bed_by_sample_to_bed(spurious_dict, bed_list, out_removed_junctions_filelist, scoring):
    sample_names = [utils.sanitize_name(bed_path) for bed_path in bed_list]
    sample_to_bed = create_sample_to_bed_dict(sample_names, out_removed_junctions_filelist)

    if out_removed_junctions_filelist is None:
        return

    sorted_keys = sorted(spurious_dict.keys())
    num_digits = len(str(len(sorted_keys)))

    for i, key in enumerate(sorted_keys):
        chrom, start, end, strand = key
        samples = spurious_dict[key]['samples']
        score2 = alignment_utils.calc_alignment_score(spurious_dict[key]['hit'], scoring)
        name = "JUNC{:0{}d}".format(i+1, num_digits)

        for _, name2, score in samples:
            sample_to_bed[name2].write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{name2}\t{score2}\n'.encode())

    for sample in sample_names:
        sample_to_bed[sample].close()
        sample_to_bed[sample] = sample_to_bed[sample].name

    return sample_to_bed


def spurious_dict_bam_by_sample_to_bed(spurious_dict, bam_list, out_removed_junctions_filelist, scoring):

    #get sample name from bam_list
    sample_names = [utils.sanitize_name(bam_path) for bam_path in bam_list]

    #dictionary where the key is the sample name and the value is a file path
    sample_to_bed = {}

    if out_removed_junctions_filelist is None:
        for sample in sample_names:
            sample_to_bed[sample] = tempfile.NamedTemporaryFile(delete=False, dir='tmp', suffix='.bed')

    else:
        sample_to_bed = create_sample_to_bed_dict(sample_names, out_removed_junctions_filelist)

    sorted_keys = sorted(spurious_dict.keys())
    named_keys = {}
    num_digits = len(str(len(sorted_keys)))
    for i, key in enumerate(sorted_keys):
        name = "JUNC{:0{}d}".format(i+1, num_digits)
        named_keys[(name,) + key] = spurious_dict[key]
        named_keys[(name,) + key]['score2'] = alignment_utils.calc_alignment_score(spurious_dict[key]['hit'], scoring)

    for (name, chrom, start, end, strand), value in named_keys.items():
        samples = value['samples']
        score2 = value['score2']
        for (sample, score) in samples:
            sample_to_bed[sample].write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{score2}\n'.encode())

    for sample in sample_names:
        sample_to_bed[sample].close()
        sample_to_bed[sample] = sample_to_bed[sample].name

    return sample_to_bed


def filter_bam_with_vacuum(bam_path, spurious_junctions_bed, out_bam_path, verbose, removed_alignments_bam):
    check_for_dependency()
    vacuum_cmd = f"{VACUUM_CMD} --remove_mate "
    if verbose:
        vacuum_cmd = f"{vacuum_cmd} -V"
    if removed_alignments_bam:
        out_bam_name, _ = os.path.splitext(out_bam_path)
        vacuum_cmd = f'{vacuum_cmd} -r {out_bam_name}_removed_alignments.bam'
    vacuum_cmd = f'{vacuum_cmd} -o {out_bam_path} {bam_path} {spurious_junctions_bed}'
    vacuum_cmd = shlex.split(vacuum_cmd)

    #use subprocess to run vacuum_cmd and return stdout and stderr
    process = subprocess.Popen((vacuum_cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if process.returncode != 0:
        print(f"Error running vacuum. Return code: {process.returncode}")
        print(f"stdout: {out.decode('utf-8')}")
        print(f"stderr: {err.decode('utf-8')}")
        sys.exit(1)

    return out.decode()

def filter_multi_bam_with_vacuum(bam_list, sample_to_bed, out_bam_list, p, verbose, removed_alignments_bam):
    #if verbose is true, make a vector of True values for each bam file
    if verbose:
        verbose = [True for bam in bam_list]
    else:
        verbose = [False for bam in bam_list]

    if removed_alignments_bam:
        removed_alignments_bam = [True for bam in bam_list]
    else:
        removed_alignments_bam = [False for bam in bam_list]


    sample_names = [utils.sanitize_name(bam_path) for bam_path in bam_list]
    #run filter_bam_with_vacuum in parallel with multiprocessing starmap
    pool = multiprocessing.Pool(processes=p)
    with pool:
        outs = pool.starmap(filter_bam_with_vacuum, zip(bam_list, [sample_to_bed[sample] for sample in sample_names],
                                                 out_bam_list, verbose, removed_alignments_bam))
    if verbose:
        for out in outs:
            print(out)

def check_for_dependency():
    """Check if runtime dependency exists."""
    if not os.path.exists(VACUUM_CMD):
        raise RuntimeError(f"{VACUUM_CMD} not found.")
