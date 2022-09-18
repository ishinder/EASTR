# Emend Alignment of Spliced Transcript Reads (EASTR)

<!-- TODO: Give a quick sentence or two to explain what this should do/give you. -->
## Required Dependencies

- [RegTools](https://regtools.readthedocs.io/en/latest/)
- [Seqtk](https://github.com/lh3/seqtk)

## Getting Started

1. Clone source repository

	```bash
	git clone https://github.com/ishinder/EASTR
	cd EASTR
	```

2. Install EASTR
	```bash
	# (OPTIONAL) Install in a Python virtual environment
	# python3 -m virtualenv venv # (OPTIONAL)
	# source ./venv/bin/activate # (OPTIONAL)
	make install # Install EASTR package
	```

## Usage
```shell
usage: EASTR [-h] [-R REFERENCE] [-bam BAM] [-A A] [-B B] [-O O O] [-E E E] [-k K] [--scoreN SCOREN] [-w W] [-p P] [-o FILE]

Emend alignments of spuriously spliced transcript reads

options:
  -h, --help            show this help message and exit
  -R REFERENCE, --reference REFERENCE
                        reference fasta genome used in alignment
  -bam BAM              Input BAM file to emend alignments
  -A A                  Matching score, default = 2
  -B B                  Mismatching penalty, default = 4
  -O O O                Gap open penalty, default = [4, 24]
  -E E E                Gap extension penalty, default = [2, 1]. A gap of length k costs min(O1+k*E1, O2+k*E2).
  -k K                  kmer length for alignment, default=7
  --scoreN SCOREN       Score of a mismatch involving ambiguous bases, default=1
  -w W                  minimizer window size, default=7
  -p P                  Number of parallel processes, default=1
  -o FILE               write output to FILE; the default output is to terminal
```

### Running a test file


<!-- TODO: Explain what output you are expected to see and brief explanation -->
```shell
eastr -R tests/data/chrX.fa \
    -bam tests/data/ERR188044_chrX.bam \
    -o tests/output/
```
