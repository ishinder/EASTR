# adapted from hisat2_extract_splice_sites.py
# https://github.com/DaehwanKimLab/hisat2/blob/master/hisat2_extract_splice_sites.py

from collections import defaultdict as dd


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
            junctions[(chrom, exons[i-1][1], exons[i][0] + 1, strand)].append(tran) #bed coordinates
    
    junctions = sorted(junctions)

    return junctions