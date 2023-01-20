import os
import tempfile
import pysam
import collections


def get_intron(alignment):
    introns = set()
    currentloc = alignment.pos
    for i,cigarop in enumerate(alignment.cigar):
        if (cigarop[0]==4): #substitution in query
            continue
        if (cigarop[0]==1): #insertion in query
            continue
        if(cigarop[0]==3):
            
            try:
                XS = alignment.get_tag('XS')
            except KeyError:
                XS = '?'
            key = (alignment.reference_name,currentloc,currentloc+cigarop[1],XS)
            introns.add(key)         
        currentloc=currentloc+cigarop[1]
    return introns

def get_spurious_alignments(bam_path, spurious_introns):
    samfile = pysam.AlignmentFile(bam_path)
    NH = collections.defaultdict((int))
    spurious_alignments = set()

    for alignment in samfile.fetch(until_eof=True):
        if not alignment.is_mapped:  # type: ignore
            continue

        introns = get_intron(alignment)
        if introns:
            for intron in introns:
                if intron in spurious_introns:
                    a = (alignment.query_name, alignment.reference_name, alignment.pos, alignment.cigarstring) 
                    if a not in spurious_alignments:
                        NH[(alignment.query_name, alignment.is_read1)] += 1
                        spurious_alignments.add(a)


    return spurious_alignments, NH


def write_filtered_bam(bam_path, outbam, spurious_introns):
    spurious_alignments, NH = get_spurious_alignments(bam_path, spurious_introns)
    removed_reads = set()
    samfile = pysam.AlignmentFile(bam_path)
    outf = pysam.AlignmentFile(outbam, "wb", template=samfile)

    for alignment in samfile.fetch(until_eof=True, multiple_iterators=True):
        if alignment.is_unmapped:
            continue 

        a = (alignment.query_name, alignment.reference_name, alignment.reference_start, alignment.cigarstring)
        if a in spurious_alignments:
            if NH[(alignment.query_name,alignment.is_read1)] == 1:
                read = 'read1' if alignment.is_read1 else 'read2'
                removed_reads.add((alignment.query_name, read))
            continue
        
        if (alignment.query_name, alignment.is_read1) in NH:
            new_NH = alignment.get_tag("NH") - NH[(alignment.query_name,alignment.is_read1)]
            alignment.set_tag("NH", new_NH)

            if new_NH <= 0:
                raise Exception("NH tag cannot be zero or negative")

        w=outf.write(alignment) 

    outf.close()
    return removed_reads



# if __name__ == '__main__':
#     # import time
#     # start = time.time()
#     # spurious_alignments, NH = get_spurious_alignments(bam_path, spurious_introns)
#     # end = time.time()
#     # print(f"took {(end-start)/60} mins"))

