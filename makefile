#!/bin/sh

all: check_make check_gcc check_input_files sureau_compiled.out 

# Check make is installed -----------------------------------------------------
check_make:
 	ifeq (, $(shell which make))
 	$(error "No make found, consider doing apt-get install make")
 	endif

# Check gcc is installed -------------------------------------------------------
check_gcc:
 	ifeq (, $(shell which gcc))
 	$(error "No gcc found, consider doing apt-get install gcc")
 	endif


# Check necessary files are available ------------------------------------------
check_input_files:
	
	@# Print working directory
	@echo Searching sureau files at:
	@pwd
	
	@# Make sure sureau.c is available
	@if [ -e ./1_source_code/sureau.c ]; then \
        	echo ""; \
        else \
        	echo "sureau.c file not found inside 1_source_code folder. Add it before continuing"; \
        	exit 1; \
	fi
	
	@# Make sure sureau_ini.txt is available	
	@if [ -e ./2_sureau_inputs/sureau_ini.txt ]; then \
        	echo ""; \
        else \
        	echo "sureau_ini.txt file not found inside 2_sureau_inputs folder. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
        	exit 2; \
	fi
	
	@# Make sure sureau_para.txt is available	
	@if [ -e ./2_sureau_inputs/sureau_para.txt ]; then \
        	echo ""; \
        else \
        	echo "sureau_para.txt file not found inside 2_sureau_inputs folder. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
        	exit 2; \
	fi
	
	@# Make sure climat_xxx_in.txt is available
	@if  [ $(shell find ./2_sureau_inputs -name 'climat*' -type f | wc -l) -eq 0 ]; then \
		echo "Climate file not found in 2_sureau_inputs. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
		echo "Climate file name should start with climat"; \
		exit 3; \
	fi
	
# Compile sureau.c -------------------------------------------------------------
#sureau.out: SHELL:=/bin/bash   # HERE: this is setting the shell for b only
sureau_compiled.out: 
	@echo Compiling sureau.c
	@gcc ./1_source_code/sureau.c -o ./2_sureau_inputs/sureau_compiled.out -lm -g
	@ echo Done! Compiled version of sureau saved at 2_sureau_inputs folder 
	@bash -c "ls"


