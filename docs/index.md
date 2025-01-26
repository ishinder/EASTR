<style>
  .md-typeset h1,
  .md-content__button {
    display: none;
  }
</style>

<p align="center">
(\(\<br>
(-.-)<br>
o('')('')<br>
</p>


<p align="center">
<em>Emending Alignment of Spliced Transcript Reads.</em>
</p>

<p align="center">
<a href="https://github.com/ishinder/EASTR/blob/main/LICENSE" target="_blank">
    <img src="https://img.shields.io/conda/l/bioconda/eastr" alt="License">
</a>
<a href="https://anaconda.org/bioconda/eastr" target="_blank">
    <img src="https://img.shields.io/conda/v/bioconda/eastr" alt="Conda Version">
</a>
<a href="https://pypi.org/project/eastr" target="_blank">
    <img src="https://img.shields.io/pypi/v/eastr" alt="PyPi">
</a>
<a href="https://anaconda.org/bioconda/eastr/files" target="_blank">
    <img src="https://img.shields.io/conda/pn/bioconda/eastr" alt="Conda Platform">
</a>
</p>


---



**EASTR** is a tool for detecting and eliminating spuriously spliced alignments in
RNA-seq datasets. It improves the accuracy of transcriptome assembly by
identifying and removing misaligned spliced alignments. The tool can process
GTF, BED, and BAM files as input. EASTR can be applied to any RNA-seq dataset
regardless of the alignment software used.

## Quickstart

Install using `conda`:

```shell
conda install bioconda:eastr
```

Install using `pip`:

```shell
pip install eastr
```

!!! warning
    Installing with pip requires you to install
    [bowtie2 >= 2.5.2](https://github.com/BenLangmead/bowtie2)
    and [samtools >= 1.19](https://github.com/samtools/samtools).


## Usage

The `tests` directory, which contains tests scripts like `run_eastr.sh` and `run_all.sh`, is available in our [GitHub repository](https://github.com/ishinder/EASTR). These scripts demonstrate two different ways to run the EASTR pipeline: using a bamlist and using a GTF file. Below, we provide detailed instructions for each use case.

### Bowtie2 Index Requirement

EASTR requires a Bowtie2 index. If you do not already have a Bowtie2 index for your reference genome, you can generate one using the following command:

```bash
bowtie2-build /path/to/reference_genome.fasta /path/to/output/bowtie2_index
```

!!! warning "Species-specific Considerations"
    - For _Homo sapiens_ and other species with pseudoautosomal regions, ensure that these regions are masked before creating a Bowtie2 index.
    - Additionally, exclude alternative scaffolds for _Homo sapiens_ to prevent inflated repeat counts.



### Running EASTR on a bamlist

1. **Ensure you are in the appropriate directory** containing the `BAM/original` folder and reference files.
2. **Create a list of BAM files** (make sure the list contains the full paths to the BAM files):

    ```bash
    ls path/to/BAM/original/*.bam > bamlist.txt
    ```

3. **Run the EASTR pipeline on the bamlist** with the following command:

    ```bash
    eastr \
        --bam bamlist.txt \
        --reference /path/to/reference_fasta \
        --bowtie2_index /path/to/bowtie2_index \
        --out_filtered_bam /path/to/output/BAM/filtered  # optional
        --out_original_junctions /path/to/output/original_junctions # optional
        --out_removed_junctions /path/to/output/removed_junctions # optional
        --removed_alignments_bam # optional
        --verbose # optional
        -p 12 # optional
    ```

### Running EASTR on a GTF

Run the EASTR pipeline on the GTF file with the following command:

```bash
  eastr
    --gtf /path/to/gtf_file
    --reference /path/to/reference_fasta
    --bowtie2_index /path/to/bowtie2_index
    --out_removed_junctions /path/to/output/outfile.bed # optional
```

## Analyzing an Example Dataset

We have included a script that demonstrates the application of the EASTR
pipeline to an _Arabidopsis_ dataset featured in our study. The
[`sra_list_arabidopsis.txt`](https://github.com/ishinder/eastr/blob/main/tests/sra_list_arabidopsis.txt) file, located in the [`tests` directory](https://github.com/ishinder/eastr/tree/main/tests) on our [GitHub repository](https://github.com/ishinder/eastr), lists the
accession IDs of the samples analyzed.

!!! warning "Required Software"
    To successfully run the example analysis, ensure that the following software is installed:

    - [`gffread >= 0.12.7`](https://github.com/gpertea/gffread)
    - [`sra-toolkit >= 3.0.1`](https://github.com/ncbi/sra-tools)

    These tools are necessary for:

    1. **Downloading FASTQ files** using the `get_fastq.sh` script.
    2. **Converting the GFF reference annotation to GTF** using the `get_ref.sh` script

The EASTR pipeline takes BAM files as input. The `run_all.sh` script acquires
FASTQ files, the FASTA reference and annotation, and then aligns the FASTQ files
using HISAT2 to generate BAM files. These BAM files are subsequently used as
input to EASTR. Additionally, EASTR can accept a GTF annotation file and output
a BED file containing questionable junctions (executed in the last command of
the `run_eastr.sh` script).

To execute the entire EASTR pipeline, which filters BAM files and identifies
reference annotation errors, use the `run_all.sh` script found in the [`tests`
directory](https://github.com/ishinder/eastr/tree/main/tests). This script ensures
all necessary steps and subscripts are carried out in the correct order. To analyze
the example dataset, follow these steps:

1. **Navigate to the `tests` directory** within the EASTR package:
2. **Make sure all scripts are executable** (`chmod +x *sh`):
3. **Run the `run_all.sh` script**.

    The script will download the necessary FASTQ files, reference genome, and then
    perform the alignment and EASTR analysis. The output files will be generated in
    their respective directories within the `tests` folder.

When executed on 4 CPUs, the EASTR command to filter 6 BAM files completes in
approximately 35 minutes, with the bulk of this time being dedicated to the
filtering of BAM files (a single BAM file typically takes between 15-20 minutes
to filter on one CPU). On 1 CPU, the EASTR command to identify questionable
introns in an annotation takes about 30 seconds.

## Command line options

```text
{{ eastr_help }}
```

## Citation

To cite EASTR in publications, please use the following reference:

!!! quote ""
    Shinder I, Hu R, Ji HJ, Chao KH, Pertea M. EASTR: Identifying and eliminating
    systematic alignment errors in multi-exon genes. Nat Commun. 2023 Nov
    9;14(1):7223. doi:
    [10.1038/s41467-023-43017-4](https://doi.org/10.1038/s41467-023-43017-4). PMID:
    [37940654](https://pubmed.ncbi.nlm.nih.gov/37940654/); PMCID:
    [PMC10632439](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10632439/).
