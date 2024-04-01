#!/bin/bash
NCPU=$1
mkdir -p ref/hisat2

REF_FA=$(ls ref/*fna)
REF_GTF=$(ls ref/*gtf)

hisat2_extract_splice_sites.py $REF_GTF > genome.ss
hisat2_extract_exons.py $REF_GTF > genome.exon
hisat2-build -p $NCPU --exon genome.exon --ss genome.ss $REF_FA ref/hisat2/hisat_index &> build_hisat_index.log &
rm genome.ss genome.exon
