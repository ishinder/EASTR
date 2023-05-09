#!/bin/bash

BASEDIR=$(pwd -P)
NCPUS=$1
OUTDIR=$BASEDIR/output/
mkdir -p $BASEDIR/BAM/filtered $BASEDIR/logs $OUTDIR

#create a bamlist of all files (input to EASTR):
find ${BASEDIR}/BAM/original/ -name "*.bam" > bamlist.txt

#bowtie2 index (NOTE: should not contain alt contigs!)
BT2_INDEX="${BASEDIR}/ref/bowtie2/bowtie2_index"
REF_FA=$(ls $BASEDIR/ref/*.fna)

# EASTR requires vacuum and junction_extractor to be in the path
# check if vacuum and junction_extractor are in the path
if ! type "vacuum" > /dev/null; then
    echo "vacuum not found in path"
    exit 1
fi

if ! type "junction_extractor" > /dev/null; then
    echo "junction_extractor not found in path"
    exit 1
fi

#run eastr on a bamlist (NOTE: bamlist should contain full paths to bam files)
(/usr/bin/time -f "\nreal\t%E\nuser\t%U\nsys\t%S" eastr \
    --verbose \
    --bam bamlist.txt \
    --reference $REF_FA \
    --bowtie2_index $BT2_INDEX \
    --out_original_junctions ${OUTDIR}/original_junctions \
    --out_removed_junctions ${OUTDIR}/removed_junctions \
    --out_filtered_bam ${BASEDIR}/BAM/filtered \
    --removed_alignments_bam \
    -p $NCPUS &> eastr_run.log) 2> time_bamlist.log


#run eastr on gtf
BASEDIR=$(pwd -P)
GTF="$(ls ${BASEDIR}/ref/*.gtf)"

(/usr/bin/time -f "\nreal\t%E\nuser\t%U\nsys\t%S" eastr \
    --gtf $GTF \
    -p $NCPUS \
    -r $REF_FA \
    -i $BT2_INDEX \
    --out_removed_junctions spurious_introns_in_gtf.bed) 2> time_gtf.log
