#!/bin/bash

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
    --option-file $sra_list_file \
    --output-directory fastq

for file in fastq/SRR*/; do
    fasterq-dump --split-files "$file" -O fastq
done

# remove subdirectories in the fastq directory
rm -r fastq/*/
