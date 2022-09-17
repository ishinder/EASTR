# Emend Alignment of Spliced Transcript Reads

## Getting Started

```bash
# Clone repository
git clone https://github.com/ishinder/EASTR
cd EASTR

# Create and activate virtual environment
python -m virtualenv .venv
source .venv/bin/activate

# Install package
pip install .
```

## Running a test file

```
eastr -R tests/data/chrX.fa \
-bam tests/data/ERR188044_chrX.bam \
-o tests/output/ERR188044_chrX_junctions.tsv 
```

## Dependencies

1. [regtools](https://regtools.readthedocs.io/en/latest/)
2. [seqtk](https://github.com/lh3/seqtk)
