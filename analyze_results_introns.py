from collections import defaultdict
from EASTR import get_spurious_introns, extract_junctions, alignment_utils
import glob
import argparse
import os

def parse_args(arglist):
    parser = argparse.ArgumentParser(
    prog="EASTR_test intron sensitivity/precision",
    description="calculate EASTR sensitivy/precision per experiment")

    parser.add_argument("--dir", help="base directory", default=os.getcwd())
    parser.add_argument("-p", help="number of processes", type=int)
    parser.add_argument("-r", help="reference GTF for sensitivity/precision calculation")

    return parser.parse_args()

def main(arglist=None):
    args = parse_args(arglist)
    p = args.p
    BASEDIR = args.dir
    ref_gtf = args.r

    bam_list = {}
    for t in ["original","filtered"]:
        bam_list[t] = glob.glob(f"{BASEDIR}/BAM/{t}/*bam")

    introns = {}
    for t in ["original","filtered"]: 
        introns[t] = get_spurious_introns.get_introns_multi_bam_regtools(bam_list[t],p)

    introns["reference"] = extract_junctions.extract_splice_sites_gtf(ref_gtf)


    for k,v in introns['filtered'].items():
        score_filt = v['score']
        score_orig = introns['original'][k]['score']

        if score_filt > score_orig:
            break

    #num of junctions in reference:
    analysis = defaultdict(dict)

    print('\n')

    for t in ["original","filtered"]:  
        num_junctions = 0

        for k,v in introns[t].items():
            num_junctions += v['score']

        analysis[t]["reference_junctions"] = set(introns[t]).intersection(set(introns['reference']))
        analysis[t]["novel_junctions"] = set(introns[t]) - set(introns['reference'])

        len_all_ref = len(introns['reference'])
        len_all = len(introns[t])
        len_x_ref = len(analysis[t]["reference_junctions"])
        len_novel = len(analysis[t]["novel_junctions"])

        sensitivity = len_x_ref/len_all_ref
        precision = len_x_ref/(len_x_ref + len_novel)


        print("****** ", t, " ******")
        print(f"total count of spliced alignments: {num_junctions}")
        print(f"reference junctions: {len_x_ref}")
        print(f"novel junctions: {len_novel}")
        print(f"sensitivity: {sensitivity:.3f}")
        print(f"precision: {precision:.3f}")
        print("\n")



    #reference junctions filtered out:
    ref_filtered = list(set(analysis['original']["reference_junctions"]) - set(analysis['filtered']["reference_junctions"]))

    for j in ref_filtered:
        prt = f"{j[0]}:{j[1]}-{j[2]} was filtered out"

        if introns['original'][j]['score'] == 1:
            prt = prt + " because it had a score of 1."
        else:
            prt = prt + "."
        print(prt)
    
    print('\n')



if __name__=="__main__":
     main()
    #settings
    # p = 12
    # BASEDIR = "/ccb/salz8-2/shinder/projects/EASTR_tests/chrX_data"
    # ref_gtf = "/ccb/salz8-2/shinder/projects/EASTR_tests/chrX_data/chrX.gtf"

