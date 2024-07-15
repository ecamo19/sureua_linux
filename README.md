# Folder Structure

The following folder structure __must__ be followed in order to run this version using the make file 

```
[root]  ./
├── [1]  1_source_code/
│   └── [311K]  sureau.c
├── [2]  2_sureau_inputs/
│   ├── [ 339]  climat_day_in.txt
│   ├── [409K]  sureau_compiled.out
│   ├── [1.2K]  sureau_ini.txt
│   └── [ 945]  sureau_para.txt
└── [3]  3_sureau_outputs/
    ├── [1005]  annual_out.csv
    ├── [4.1K]  sureau_out.csv
    └── [8.1M]  transient_out.csv
```
# Running SurEau

For running sureau just type `make run_sureau` in the terminal console. GCC and make must be installed prior to running the model

# Example run

```bash
cd example_run
make example_run
```
