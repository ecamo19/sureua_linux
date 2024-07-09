
all: check_make check_gcc check_input_files sureau.out 

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
        	echo "File found"; \
        else \
        	echo "File sureau.c not found inside 1_source_code folder. Add it before continuing"; \
        	exit; \
	fi
# Compile sureau.c -------------------------------------------------------------

# if folder 1_source_code and sureau.c not found
# Create folder and download sureau.c from repo

# if folder 1_source_code not found but sureau.c found
# Create 1_source_code and move sureau.c inside

# else
# compile code

#sureau.out: SHELL:=/bin/bash   # HERE: this is setting the shell for b only

sureau_compiled.out: 
	@echo Compiling Sureau
	@gcc -c ./1_source_code/sureau.c -o ./2_sureau_inputs/sureau_compiled.out -lm -g
	@bash -c "ls"
# create folder struture

