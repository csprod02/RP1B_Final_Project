# RP1B_Final_Project

Included in this repository are the three Python scripts that were used to generate the results discussed in this project.

## rp1b_script.py

`rp1b_script.py` was initially run on the data included in the `refSeqs_and_VCF` directory.

The script takes the following arguments as input from the command line:
- `-i` is the input multisample VCF file path
- `-r` is the input reference genome file path (FASTA format for example)
- `-j` is the job title which is used to distinguish files generated when running on multiple occasions
- `-w` can be optionally used to alter the window size (measured in Kb) used during the identification of recombinant SNPs (the default value is 2, representing a window size of 2000 bases)
- `-s` can be optionally used to alter the threshold of SNPs required in the given window when deciding whether to classify the region as recombinant or not

## rp1b_plots.py

`rp1b_plots.py` plots the graphs shown as Figures 1-5 of my written report.

This script was designed to be run specifically on the data produced during this project from `rp1b_script.py`.  Therefore, no arguments are required.  **However**, the job title used during the running of `rp1b_script.py` must be consistent with the file names used here.  **Therefore** it is suggested that when running `rp1b_script.py` on the data again, the same job titles are used and that the commands follow what is advised previously.

## rp1b_dim_reduction.py
