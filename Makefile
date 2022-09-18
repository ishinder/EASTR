PACAKGE = EASTR
.PHONY: all clean

all: install

install:
	@echo "Checking if RegTools is installed..."
	@which regtools || echo "RegTools not found. Please install";exit 1
	@echo "Checking if Seqtk is installed..."
	@which seqtk || echo "Seqtk not found. Please install"; exit 1
	@pip install .

clean:
	@echo "Cleaning"
	@pip uninstall -y ${PACAKGE}
	@echo "Done"
