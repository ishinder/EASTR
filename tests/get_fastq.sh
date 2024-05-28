#!/bin/bash
set -xe

# create fastq directory
mkdir -p fastq

while getopts f: flag; do
    case "${flag}" in
        f) sra_list_file=${OPTARG};;
        *) echo "Invalid";;
    esac
done

# download and extract data into fastq directory
prefetch -p \
    --option-file "$sra_list_file" \
    --output-directory fastq

for file in fastq/SRR*/; do
    fasterq-dump -p -v --split-files "$file" --outdir fastq
done

# remove subdirectories in the fastq directory
rm -r fastq/*/
