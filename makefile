run_sureau: check_make check_gcc check_input_files sureau_compiled.out 

# Check make is installed -----------------------------------------------------
check_make:
 	ifeq (, $(shell which make))
 	$(error "No make found, consider doing apt-get install make")
 	endif

# Check gcc is installed -------------------------------------------------------
check_gcc:
 	ifeq (, $(shell which gcc))
 	$(error "N
 	o gcc found, consider doing apt-get install gcc")
 	endif


# Check necessary files are available ------------------------------------------
check_input_files:
	
	@# Print working directory
	@echo Searching sureau files at:
	@pwd
	
	@# Make sure sureau.c is available
	@if [ -e ./1_source_code/sureau.c ]; then \
        	echo "\nsureau_c_file found\n"; \
        else \
        	echo "sureau.c file not found inside 1_source_code folder. Add it before continuing"; \
        	exit 1; \
	fi
	
	@# Make sure sureau_ini.txt is available	
	@if [ -e ./2_sureau_inputs/sureau_ini.txt ]; then \
        	echo "sureau_init.txt found\n"; \
        else \
        	echo "sureau_ini.txt file not found inside 2_sureau_inputs folder. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
        	exit 2; \
	fi
	
	@# Make sure sureau_para.txt is available	
	@if [ -e ./2_sureau_inputs/sureau_para.txt ]; then \
        	echo "sureau_para.txt file found\n"; \
        else \
        	echo "sureau_para.txt file not found inside 2_sureau_inputs folder. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
        	exit 2; \
	fi
	
	@# Make sure climat_xxx_in.txt is available
	@if  [ $(shell find ./2_sureau_inputs -name 'climat*' -type f | wc -l) -eq 0 ]; then \
		echo "Climate file not found in 2_sureau_inputs. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
		echo "Climate file name should start with climat"; \
		exit 3; \
	elif [ $(shell find ./2_sureau_inputs -name 'climat*' -type f | wc -l) -eq 1 ]; then \
		echo "climat file found\n"; \
	elif [ $(shell find ./2_sureau_inputs -name 'climat*' -type f | wc -l) -ge 2 ]; then \
		echo "Two or more climate files found. Use just one"; \
		exit 4; \
	else \
		echo "Failed reading climat files"; \
		exit 5; \
	fi
	
# Compile sureau.c -------------------------------------------------------------
#sureau.out: SHELL:=/bin/bash   # HERE: this is setting the shell for b only
sureau_compiled.out: 

	@gcc ./1_source_code/sureau.c -o ./2_sureau_inputs/sureau_compiled.out -lm -g
	@ echo Done! Compiled version of sureau saved at 2_sureau_inputs folder 

# Run SurEau -------------------------------------------------------------------
	@if [ -e ./2_sureau_inputs/sureau_out.csv ]; then \
        	echo "\nsureau_out.csv found\n"; \
        else \
        	echo "sureau_out.csv file not found inside 2_sureau_inputs folder."; \
			touch ./2_sureau_inputs/sureau_out.csv; \
			cd 2_sureau_inputs; \
			pwd; \
	fi
	
# If sureau_out.csv does not exist create a new empty sureau_out.csv with touch

