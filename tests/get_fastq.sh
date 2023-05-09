#!/bin/bash

# create fastq directory
mkdir fastq

while getopts f: flag; do
    case "${flag}" in
        f) sra_list_file=${OPTARG};;
    esac
done

# download and extract data into fastq directory
prefetch --option-file $sra_list_file --output-directory fastq

for i in $(ls -d fastq/SRR*/); do 
    fasterq-dump --split-files $i -O fastq
done

# remove subdirectories in the fastq directory
rm -r fastq/*/