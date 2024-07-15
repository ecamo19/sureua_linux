# Run all rules -----------------------------------------------------------------
run_sureau: check_1_source_code_folder check_2_sureau_inputs_folder \
			check_3_sureau_outputs_folder \
			check_make check_gcc check_input_files sureau_compiled.out 

# Check make is installed -------------------------------------------------------
check_make:
 	ifeq (, $(shell which make))
 	$(error "No make found, consider doing apt-get install make")
 	endif

# Check gcc is installed --------------------------------------------------------
check_gcc:
 	ifeq (, $(shell which gcc))
 	$(error "N
 	o gcc found, consider doing apt-get install gcc")
 	endif

# Chech folders exists ----------------------------------------------------------
check_1_source_code_folder:
	@if [ -d ./1_source_code ]; then\
  		echo "Directory 1_source_code found";\
	else \
		echo "\nCreate folder 1_source_code before continuing\n"; \
		exit 1;\
	fi

check_2_sureau_inputs_folder:
	@if [ -d ./2_sureau_inputs ]; then\
  		echo "Directory 2_sureau_inputs found";\
	else \
		echo "\nCreate folder 2_sureau_inputs before continuing\n"; \
		exit 2;\
	fi

check_3_sureau_outputs_folder:
	@if [ -d ./3_sureau_outputs ]; then\
  		echo "Directory 3_sureau_outputs found";\
	else \
		echo "\nCreate folder 3_sureau_outputs before continuing\n"; \
		exit 3;\
	fi

# Check the necessary files are available ---------------------------------------
check_input_files:
	
	@# Print working directory
	@echo Searching sureau files at:
	@pwd
	
	@# Make sure sureau.c is available
	@if [ -e ./1_source_code/sureau.c ]; then \
        	echo "\nsureau_c_file found\n"; \
    else \
        	echo "sureau.c file not found inside 1_source_code folder. Add it before continuing"; \
        	exit 4; \
	fi
	
	@# Make sure sureau_ini.txt is available	
	@if [ -e ./2_sureau_inputs/sureau_ini.txt ]; then \
        	echo "sureau_init.txt found\n"; \
        else \
        	echo "sureau_ini.txt file not found inside 2_sureau_inputs folder. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
        	exit 5; \
	fi
	
	@# Make sure sureau_para.txt is available	
	@if [ -e ./2_sureau_inputs/sureau_para.txt ]; then \
        	echo "sureau_para.txt file found\n"; \
        else \
        	echo "sureau_para.txt file not found inside 2_sureau_inputs folder. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
        	exit 6; \
	fi
	
	@# Make sure climat_xxx_in.txt is available
	@if  [ $(shell find ./2_sureau_inputs -name 'climat*' -type f | wc -l) -eq 0 ]; then \
		echo "Climate file not found in 2_sureau_inputs. Use the sureau_parameters_cheatsheet.xlsx file for creating one"; \
		echo "Climate file name should start with climat"; \
		exit 7; \
	elif [ $(shell find ./2_sureau_inputs -name 'climat*' -type f | wc -l) -eq 1 ]; then \
		echo "climat file found\n"; \
	elif [ $(shell find ./2_sureau_inputs -name 'climat*' -type f | wc -l) -ge 2 ]; then \
		echo "Two or more climate files found. Use just one"; \
		exit 8; \
	else \
		echo "Failed reading climat files"; \
		exit 9; \
	fi
	
# Compile sureau.c --------------------------------------------------------------
sureau_compiled.out: 

	@gcc ./1_source_code/sureau.c -o ./2_sureau_inputs/sureau_compiled.out -lm -g
	@echo Done! Compiled version of sureau saved at 2_sureau_inputs folder 

# Run SurEau --------------------------------------------------------------------
	@if [ -e ./2_sureau_inputs/sureau_out.csv ]; then \
        	echo "\nsureau_out.csv found\n"; \
			cd 2_sureau_inputs; \
			./sureau_compiled.out; \
			mv -t ../3_sureau_outputs annual_out.csv sureau_out.csv transient_out.csv;\
        else \
        	echo "sureau_out.csv file not found inside 2_sureau_inputs folder."; \
			echo "creating empty sureau_out.csv"; \
			touch ./2_sureau_inputs/sureau_out.csv; \
			cd 2_sureau_inputs; \
			./sureau_compiled.out; \
			mv -t ../3_sureau_outputs annual_out.csv sureau_out.csv transient_out.csv;\
	fi