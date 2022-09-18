.PHONY: all regtools seqtk clean

all: regtools seqtk

check_regtools:
	@echo "Checking if regtools is installed..."
	@which regtools

check_seqtk:
	@echo "Checking if seqtk is installed..."
	@which seqtk

install:
	@pip install .

regtools:
	@echo "Building regtools"
	@cd regtools/; mkdir build; cd build/; cmake ..; make

seqtk:
	@echo "Building seqtk"
	@$(MAKE) -C seqtk

clean:
	@echo "Cleaning"
	@$(MAKE) -C regtools/build clean
	@$(MAKE) -C seqtk/ clean
	@echo "Done"
