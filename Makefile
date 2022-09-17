.PHONY: all regtools seqtk clean

all: regtools seqtk


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
