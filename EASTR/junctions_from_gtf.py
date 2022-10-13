#adapted from hisat2_extract_splice_sites.py

from collections import defaultdict as dd
import pandas as pd
import glob

def extract_splice_sites(gtf_path):
    genes = dd(list)
    trans = {}

    gtf =open(gtf_path, "r")

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf:
        line = line.strip()
        if line.startswith('#'):
            continue
        chrom, source, feature, start, end, score, \
                strand, frame, attributes = line.split('\t')

        start, end = int(start), int(end)

        if start > end:
            raise Exception("Start of region can not be greater than end of region for:\n",line)

        if feature != 'exon':
            continue

        values_dict = {}
        for attr in attributes.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            raise Exception("Exon does not contain transcript or gene ID\n")

        transcript_id = values_dict['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[start, end]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            trans[transcript_id][2].append([start, end])


    for tran, [chrom, strand, exons] in trans.items():
            exons.sort()


    junctions = dd(list)
    for tran, (chrom, strand, exons) in trans.items():
        for i in range(1, len(exons)):
            junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)].append(tran) #intron bed coordinates
    
    return junctions

# ref_gtf = '/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf'
# ref_junctions = extract_splice_sites(ref_gtf)
# cols = ['seqname','start','end','strand','transcripts']
# ref_junctions = pd.Series(ref_junctions).rename_axis(cols[0:4]).reset_index(name=cols[4])


# files = {"m40w7k7":[],"m14w1k3":[]}
# for setting in files:

#     f = glob.glob(f"/ccb/salz8-2/shinder/projects/EASTR_tests/chess_brain/BED/R2816*{setting}*.bed")
#     files[setting] = f

# spurious = {"m40w7k7":pd.DataFrame(),"m14w1k3":pd.DataFrame()}

# spurious["m40w7k7"]=pd.read_csv(files["m40w7k7"][0],header=None,sep='\t', names=['seqname','start','end',"AS"])
# spurious["m14w1k3"]=pd.read_csv(files["m14w1k3"][0],header=None,sep='\t', names=['seqname','start','end',"AS"])

# compare = {}
# compare["m40w7k7"] = pd.merge(spurious["m40w7k7"],ref_junctions,how="inner",on=['seqname','start','end'])
# compare["m14w1k3"] = pd.merge(spurious["m14w1k3"],ref_junctions,how="inner",on=['seqname','start','end'])

