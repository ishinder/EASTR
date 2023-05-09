
# Emending Alignments of Spliced Transcript Reads (EASTR)
\(\\(\\  
\(-.-\)  
o\(''\)(''\)  

EASTR is a tool for detecting and removing spurious splice junctions in RNA-seq datasets. It improves the accuracy of transcriptome assembly and quantification by identifying and removing misaligned spliced alignment. The tool can process GTF, BED, and BAM files as input. EASTR can be applied to any RNA-seq dataset regardless of the alignment software used.


## Required Dependencies
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)
- [mappy: Minimap2 Python Binding](https://github.com/lh3/minimap2/tree/master/python) 

## Getting Started
The installation steps for running EASTR are outlined below. Installation takes a few minutes. 

1. Clone source repository
	```bash
	git clone --recursive https://github.com/ishinder/EASTR.git
	cd EASTR
	```
2. Compile junction_extractor and vacuum
   (Note: )
    ```bash
    cd utils
    cmake
    make
    ```

3. To add the junction_extractor and vacuum executables in the utils folder to your PATH, you can follow these steps:  

    Get the absolute path of the build directory:
    ```bash 
    pwd
    ```
    This command will print the absolute path of the current directory, which is the `utils` directory.

    Add the absolute path to the PATH environment variable:  
    For temporary use (not persistent across sessions), you can run:

    ```bash
    export PATH=$PATH:<absolute_path_to_utils_directory>
    ```

    Replace `<absolute_path_to_utils_directory>` with the path you obtained in the previous step.
 
    For persistent use (across sessions), you can add the `export` command to your shell's configuration file. For bash, this is typically the `.bashrc` or `.bash_profile` file in your home directory. For zsh, this is the `.zshrc` file. To add the path to your shell configuration file, run:

    ```bash
    echo 'export PATH=$PATH:<absolute_path_to_utils_directory>' >> ~/.bashrc
    ```

4. Install EASTR
	```bash
	# (OPTIONAL) Install in a Python virtual environment
	# python3 -m virtualenv .venv # (OPTIONAL)
	# source .venv/bin/activate # (OPTIONAL)
	make install # Install EASTR package
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

The `run_eastr.sh` script in the `tests` directory demonstrates two different ways to run the EASTR pipeline: on a bamlist and on a GTF file. Below, we provide instructions for each use case.

### Running EASTR on a bamlist
1. Ensure you are in the appropriate directory containing the BAM/original folder and reference files.
2. Create a list of BAM files (make sure the list contains the full paths to the BAM files):
    ```bash
    ls path/to/BAM/original/*.bam > bamlist.txt
    ```
3. Run the EASTR pipeline on the bamlist with the following command:
    ```bash
    eastr 
        --bam bamlist.txt 
        --reference /path/to/reference_fasta 
        --bowtie2_index /path/to/bowtie2_index 
        --out_filtered_bam /path/to/output/BAM/filtered  #optional
        --out_original_junctions /path/to/output/original_junctions #optional
        --out_removed_junctions /path/to/output/removed_junctions # optional 
        --removed_alignments_bam #optional
        --verbose #optional
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


## Analyzing an example dataset
We have included a script that demonstrates the application of the EASTR pipeline to an *Arabidopsis* dataset featured in our study. The `sra_list_arabidopsis.txt` file, located in the `tests` directory, lists the accession IDs of the samples analyzed.

The EASTR pipeline takes BAM files as input. The `run_all.sh` script acquires FASTQ files, the FASTA reference and annotation, and then aligns the FASTQ files using HISAT2 to generate BAM files. These BAM files are subsequently used as input to EASTR. Additionally, EASTR can accept a GTF annotation file and output a BED file containing questionable junctions (executed in the last command of the `run_eastr.sh` script).

To execute the entire EASTR pipeline, which filters BAM files and identifies reference annotation errors, use the `run_all.sh` script found in the `tests` directory. This script ensures all necessary steps and subscripts are carried out in the correct order. To analyze the example dataset, follow these steps:

1. Navigate to the `tests` directory within the EASTR package:
2. Make sure all scripts are executable (`chmod +x *sh`):
3. Run the `run_all.sh` script.

The script will download the necessary FASTQ files, reference genome, and then perform the alignment and EASTR analysis. The output files will be generated in their respective directories within the `tests` folder.

When executed on 4 CPUs, the EASTR command to filter 6 BAM files completes in approximately 35 minutes, with the bulk of this time being dedicated to the filtering of BAM files \(a single bam file typically takes between 15-20 minutes to filter on a single CPU). On 1 CPU, the EASTR command to identify questionable introns in an annotation takes about 30 seconds.

*Note 1*: Downloading FASTQ files using the `get_fastq.sh` script requires [SRA_toolkit](https://github.com/ncbi/sra-tools)  
*Note 2*: Converting the GFF reference annotation to GTF in the `get_ref.sh` script required [gffread](https://github.com/gpertea/gffread)

<!-- 
# Citation
To cite EASTR in publications, please use:
[citation] -->
