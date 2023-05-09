#!/bin/bash

# Get the NCBI accession number from the command line
accession_number=$1

# Download the reference from NCBI
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${accession_number}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF&filename=${accession_number}.zip" -H "Accept: application/zip"
unzip ${accession_number}.zip
rm ${accession_number}.zip

mkdir ref
rm ncbi_dataset/data/assembly_data_report.jsonl
rm ncbi_dataset/data/dataset_catalog.json
rm README.md
mv ncbi_dataset/data/${accession_number}/* ref
awk '{ if ($7 == "?") $7 = "."; print }' ref/genomic.gff > ref/genomic_strand_fixed.gff
gffread ref/genomic_strand_fixed.gff -T -o ref/genomic.gtf
rm ref/genomic_strand_fixed.gff ref/genomic.gff
