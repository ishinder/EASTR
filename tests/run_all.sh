#!/bin/bash

# Input arguments
sra_list_file="sra_list_arabidopsis.txt"
accession_number="GCA_000001735.2"
NCPU=16

# 1. get_fastq.sh
bash get_fastq.sh -f $sra_list_file

# 2. get_ref.sh
bash get_ref.sh $accession_number

# 3. build_hisat_index
bash build_hisat_index.sh $NCPU

# 4. build_bowtie_index.sh
bash build_bowtie_index.sh $NCPU

# 5. align with hisat2
bash align_hisat2.sh $NCPU

# 6. run eastr
bash run_eastr.sh $NCPU
