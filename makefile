
all: check_make check_gcc check_input_files sureau_compiled.out 

# Check if make is installed -------------------------------------------------
check_make:
 	ifeq (, $(shell which make))
 	$(error "No make found, consider doing apt-get install make")
 	endif


check_gcc:
 	ifeq (, $(shell which gcc))
 	$(error "No gcc found, consider doing apt-get install gcc")
 	endif


# First, check necessary files are available ---------------------------------
check_input_files:

	@# Print working directory
	@echo Searching the sureau_inputs folder at:
	@pwd
	
	@# Check that 2_sureau_inputs folder exist
	@if [ -e ./1_source_code/sureau.c ]; then \
        	echo ""; \
        else \
        	echo "sureau.c file not found inside 1_source_code folder. Add it before continuing"; \
        	exit; \
	fi
# Compile sureau.c -------------------------------------------------------------
#sureau.out: SHELL:=/bin/bash   # HERE: this is setting the shell for b only
sureau_compiled.out: 
	@echo Compiling sureau.c
	@gcc -c ./1_source_code/sureau.c -o ./2_sureau_inputs/sureau_compiled.out -lm -g
	@ echo Done! Compiled version of sureau saved at 2_sureau_inputs folder 
	@bash -c "ls"


