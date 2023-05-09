#!/bin/bash
NCPU=$1
mkdir -p ref/bowtie2
REF_FA=$(ls ref/*fna)

bowtie2-build --threads $NCPU $REF_FA ref/bowtie2/bowtie2_index &> build_bowtie_index.log &
