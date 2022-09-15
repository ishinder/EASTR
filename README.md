# Emend Alignment of Spliced Transcript Reads

## <a name="installation"></a>Download and Installation
```sh
#download
git clone https://github.com/ishinder/EASTR
cd EASTR

#create and activate virtual environment
virtualenv venv
. env/bin/activate

#install package
pip install .
```
## <a name="testing"></a>Running a test file
```
eastr -R tests/data/chrX.fa -bam tests/data/ERR188044_chrX.bam -o tests/output/
```

## <a name="dependencies"></a> Software Dependencies
1. [regtools](https://regtools.readthedocs.io/en/latest/commands/junctions-extract/)
2. [seqtk](https://github.com/lh3/seqtk)