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

The uvicorn command line tool is the easiest way to run your application.

### Command line options

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
