PACAKGE = EASTR
.PHONY: all clean

all: install

install:
	@echo "Checking if bowtie2 is installed..."
	@which bowtie2 || echo "bowtie2 not found. Please install"; exit 0
	@echo "Checking if samtools is installed..."
	@which samtools || echo "samtools not found. Please install"; exit 0
	@echo "Checking if junction_extractor is installed..."
	@which junction_extractor || echo "junction_extractor not found. Please install"; exit 0
	@echo "Checking if vacuum is installed..."
	@which vacuum || echo "vacuum not found. Please install"; exit 0
	@pip install .

clean:
	@echo "Cleaning"
	@pip uninstall -y ${PACAKGE}
	@echo "Done"
