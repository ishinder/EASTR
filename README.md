
# Emending Alignments of Spliced Transcript Reads (EASTR)
\(\\(\\  
\(-.-\)  
o\(''\)(''\)  

EASTR is a tool for detecting and removing spurious splice junctions in RNA-seq datasets. It improves the accuracy of transcriptome assembly and quantification by correcting misaligned reads. The tool can process GTF, BED, and BAM files as input. EASTR can be applied to any RNA-seq dataset regardless of the alignment software used.

<!-- TODO: Give a quick sentence or two to explain what this should do/give you. -->
## Required Dependencies

- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)
- [mappy: Minimap2 Python Binding](https://github.com/lh3/minimap2/tree/master/python) 

## Getting Started

1. Clone source repository

	```bash
	git clone --recursive https://github.com/ishinder/EASTR.git
	cd EASTR
	```
2. Compile junction_extractor and vacuum

    ```bash
    cd utils
    mkdir build
    cd build
    cmake ..
    make
    ```

3. To add the junction_extractor and vacuum executables in the utils folder to your PATH, you can follow these steps:  

    Get the absolute path of the build directory:
    ```bash 
    pwd
    ```
    This command will print the absolute path of the current directory, which is the `build` directory inside the `utils` folder.

    Add the absolute path to the PATH environment variable:  
    For temporary use (not persistent across sessions), you can run:

    ```bash
    export PATH=$PATH:<absolute_path_to_build_directory>
    ```

    Replace `<absolute_path_to_build_directory>` with the path you obtained in the previous step.
 
    For persistent use (across sessions), you can add the `export` command to your shell's configuration file. For bash, this is typically the `.bashrc` or `.bash_profile` file in your home directory. For zsh, this is the `.zshrc` file. To add the path to your shell configuration file, run:

    ```bash
    echo 'export PATH=$PATH:<absolute_path_to_build_directory>' >> ~/.bashrc
    ```

4. Install EASTR
	```bash
	# (OPTIONAL) Install in a Python virtual environment
	# python3 -m virtualenv venv # (OPTIONAL)
	# source .venv/bin/activate # (OPTIONAL)
	make install # Install EASTR package
	```

5. OPTIONAL: Download test data
  ```bash
  wget ftp://ftp.ccb.jhu.edu/pub/ishinder/chrX.gtf
  ```

### Required Arguments

Note: Only one of the above input options (GTF, BED, or BAM) should be provided.  
- `--gtf` : Input GTF file containing transcript annotations
- `--bed` : Input BED file with intron coordinates
- `--bam` : Input BAM file or a TXT file containing a list of BAM files with read alignments

  
Additionally, the following arguments are required:
- `-r`, `--reference` : Reference FASTA genome used in alignment
- `-i`, `--bowtie2_index` : Path to Bowtie2 index for the reference genome

### Optional Arguments

- `--bt2_k` : Minimum number of distinct alignments found by bowtie2 for a junction to be considered spurious. Default: 10
- `-o` : Length of the overhang on either side of the splice junction. Default: 50
- `-a` : Minimum required anchor length in each of the two exons. Default: 7
- `--min_junc_score` : Minimum number of supporting spliced reads required per junction. Default: 1
- `--trusted_bed` : Path to a BED file path with trusted junctions, which will not be removed by EASTR.
- `--verbose` : Display additional information during BAM filtering, including the count of total spliced alignments and removed alignments
- `--removed_alignments_bam` : Write removed alignments to a BAM file
- `-p` : Number of parallel processes. Default: 1

### Minimap2 Parameters

- `-A` : Matching score. Default: 3
- `-B` : Mismatching penalty. Default: 4
- `-O` : Gap open penalty. Default: [12, 32]
- `-E` : Gap extension penalty. Default: [2, 1]
- `-k` : K-mer length for alignment. Default: 3
- `--scoreN` : Score of a mismatch involving ambiguous bases. Default: 1
- `-w` : Minimizer window size. Default: 2
- `-m` : Discard chains with chaining score. Default: 25

### Output Options

- `--out_original_junctions` : Write original junctions to the output file or directory
- `--out_removed_junctions` : Write removed junctions to the output file or directory; the default output is to the terminal
- `--out_filtered_bam` : Write filtered bams to the output file or directory
- `--filtered_bam_suffix` : Suffix added to the name of the output BAM files. Default: '_EASTR_filtered'




## Usage



## Running a test file


<!-- TODO: Explain what output you are expected to see and brief explanation -->
```shell
eastr -R tests/data/chrX.fa \
    -bam tests/data/ERR188044_chrX.bam \
    --out_introns tests/output/ERR188044_chrX_spurious_introns.bed
    --out_bam tests/output/ERR188044_chrX_filtered.bam
    --removed_reads tests/output/ERR188044_chrX_removed_reads.txt
```

# Citation
To cite EASTR in publications, please use:
[citation]
