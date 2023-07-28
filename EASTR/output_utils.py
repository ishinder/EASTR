import csv
import shlex
import subprocess
from typing import List, Union
from io import StringIO
import multiprocessing
import os
import tempfile
from EASTR import utils, alignment_utils
import sys

def out_junctions_filelist(bam_list:list, gtf_path, bed_list, out_junctions, suffix="") -> Union[List[str], None, str]:
    if out_junctions is None:
        return None
    
    if gtf_path:
        if utils.check_directory_or_file(out_junctions) == 'dir':
            out_junctions= out_junctions + "/" + os.path.splitext(os.path.basename(gtf_path)[0]) + suffix + ".bed"     
        return out_junctions
    
    if bed_list:
        if suffix in ["_original_junctions", ""]:
            return None
        
        if len(bed_list) == 1:
            if utils.check_directory_or_file(out_junctions) == 'dir':
                out_junctions = f"{out_junctions}/{os.path.splitext(os.path.basename(gtf_path))[0]}{suffix}.bed"
            path = os.path.dirname(out_junctions)
            utils.make_dir(path)
            return [out_junctions]

        if utils.check_directory_or_file(out_junctions) == 'file':
            print("ERROR: the path provided for the output bed files is a file path, not a directory")
            sys.exit(1)
        
        utils.make_dir(out_junctions)
        out_junctions_filelist = []
        for bed in bed_list:
            out_junctions_filelist.append(f"{out_junctions}/{os.path.splitext(os.path.basename(bed))[0]}{suffix}.bed")
        return out_junctions_filelist
    
    if len(bam_list) == 1:
        if utils.check_directory_or_file(out_junctions) == 'dir':
            out_junctions = f"{out_junctions}/{os.path.splitext(os.path.basename(bam_list[0]))[0]}{suffix}.bed"
        path = os.path.dirname(out_junctions)
        utils.make_dir(path)
        return [out_junctions]

    if utils.check_directory_or_file(out_junctions) == 'file':
        print("ERROR: the path provided for the output bed files is a file path, not a directory")
        sys.exit(1)

    utils.make_dir(out_junctions)
    out_junctions_filelist = []
    for bam in bam_list:
        out_junctions_filelist.append(out_junctions + "/" + os.path.splitext(os.path.basename(bam))[0] + suffix + ".bed")
    return out_junctions_filelist

def out_filtered_bam_filelist(bam_list:list, out_filtered_bam, suffix="_EASTR_filtered") -> Union[List[str], None]:
    if bam_list is None or out_filtered_bam is None:
        return

    if len(bam_list) == 1:
        if utils.check_directory_or_file(out_filtered_bam) == 'dir':
            out_filtered_bam = out_filtered_bam + "/" + os.path.splitext(os.path.basename(bam_list[0]))[0] + suffix + ".bam"
        path = os.path.dirname(out_filtered_bam)
        utils.make_dir(path)
        out_filtered_bam_filelist = [out_filtered_bam]

    else:
        if utils.check_directory_or_file(out_filtered_bam) == 'file':
            print("ERROR: the path provided for the output file is a file, not a directory")
            sys.exit(1)
        utils.make_dir(out_filtered_bam)
        out_filtered_bam_filelist = []
        for bam in bam_list:
            out_filtered_bam_filelist.append(out_filtered_bam + "/" + os.path.splitext(os.path.basename(bam))[0] + suffix + ".bam")

    return out_filtered_bam_filelist

def writer_spurious_dict_bam_to_bed(spurious_dict, named_keys, scoring, writer):
    for key, value in spurious_dict.items():
        chrom , start, end, strand = key
        name = named_keys[key]
        score = value['score']
        samples = value['samples']
        score2 = alignment_utils.calc_alignment_score(value['hit'],scoring)
        name2 = ';'.join([f"{score2},{sample_id}" for score2, sample_id in samples])
        writer.writerow([chrom, start, end, name, score, strand, score2, name2])


def writer_spurious_dict_gtf_to_bed(spurious_dict, named_keys, scoring, writer):
    for key, value in spurious_dict.items():
        chrom , start, end, strand = key
        gene_id = value['transcripts'][0]
        tanscripts = value['transcripts'][1]
        name = named_keys[key]
        score = '.'
        score2 = alignment_utils.calc_alignment_score(value['hit'],scoring)
        name2 = ';'.join([f"{gene_id}", *tanscripts])
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
    out = StringIO()
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


def spurious_dict_bed_by_sample_to_bed(spurious_dict, bed_list, out_removed_junctions_filelist, scoring):
    sample_names = [os.path.splitext(os.path.basename(bed_path))[0] for bed_path in bed_list]
    sample_to_bed = {}
    
    if out_removed_junctions_filelist is None:
        return

    else:
        for sample in sample_names:
            for i, sample_id in enumerate(out_removed_junctions_filelist):
                if sample_id == sample:
                    file_path = out_removed_junctions_filelist[i]
                    file_obj = open(file_path, mode='w+b')
                    sample_to_bed[sample] = file_obj
                    break

    sorted_keys = sorted(spurious_dict.keys())
    for key in sorted_keys:
        chrom, start, end, strand = key
        samples = spurious_dict[key]['samples']
        score2 = alignment_utils.calc_alignment_score(spurious_dict[key]['hit'], scoring)

        for name, _, score in samples:
            sample_to_bed[name].write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{score2}\n'.encode())

    for sample in sample_names:
        sample_to_bed[sample].close()
        sample_to_bed[sample] = sample_to_bed[sample].name

    return sample_to_bed


def spurious_dict_bam_by_sample_to_bed(spurious_dict, bam_list, out_removed_junctions_filelist, scoring):

    #get sample name from bam_list
    sample_names = [os.path.splitext(os.path.basename(bam_path))[0] for bam_path in bam_list]

    #dictionary where the key is the sample name and the value is a file path
    sample_to_bed = {}
    
    if out_removed_junctions_filelist is None:
        for sample in sample_names:
            sample_to_bed[sample] = tempfile.NamedTemporaryFile(delete=False, dir='tmp', suffix='.bed')

    else:
        for sample in sample_names:
            longest_match = None
            for i, sample_id in enumerate(out_removed_junctions_filelist):
                if sample in sample_id:
                    if longest_match is None or len(sample_id) > len(longest_match):
                        longest_match = sample_id
            if longest_match is not None:  # Only update if a match was found
                file_path = out_removed_junctions_filelist[out_removed_junctions_filelist.index(longest_match)]
                file_obj = open(file_path, mode='w+b')
                sample_to_bed[sample] = file_obj

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
    vacuum_cmd = "vacuum --remove_mate "
    if verbose:
        vacuum_cmd = f"{vacuum_cmd} -V"
    if removed_alignments_bam:
        out_bam_name = os.path.splitext(out_bam_path)[0]
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


    sample_names = [os.path.splitext(os.path.basename(bam_path))[0] for bam_path in bam_list]
    #run filter_bam_with_vacuum in parallel with multiprocessing starmap
    pool = multiprocessing.Pool(processes=p)
    with pool:
        outs = pool.starmap(filter_bam_with_vacuum, zip(bam_list, [sample_to_bed[sample] for sample in sample_names], 
                                                 out_bam_list, verbose, removed_alignments_bam))
    if verbose:
        for out in outs:
            print(out)

#def write_gtf_to_bed(spurious_dict, out_removed_junctions_filelist, scoring):
    # sorted_keys = spurious_dict.keys()
    # named_keys = {}
    # for i, key in enumerate(sorted_keys):
    #     name = "JUNC{}".format(i+1)
    #     named_keys[key] = name

    # out_io_dict = {}
    # for i,sample in enumerate(sample_names):
    #     out_io_dict[sample] = [StringIO(), out_removed_junctions_filelist[i]]

    # for key, value in spurious_dict.items():
    #     name = named_keys[key]
    #     samples = value['samples']
    #     score2 = alignment_utils.calc_alignment_score(value['hit'],scoring)
    #     for i, sample in enumerate(samples):
    #         score, sample_id = sample
    #         out_io_dict[sample_id][0].write(f"{key[0]}\t{key[1]}\t{key[2]}\t{name}\t{score}\t{key[3]}\t{score2}")

    # for (out,filepath) in out_io_dict.values():
    #     with open(filepath, 'w') as out_removed_junctions:
    #         out_removed_junctions.write(out.getvalue())

# if __name__ == '__main__':
#     # import time
#     # start = time.time()
#     # spurious_alignments, NH = get_spurious_alignments(bam_path, spurious_introns)
#     # end = time.time()
#     # print(f"took {(end-start)/60} mins"))

