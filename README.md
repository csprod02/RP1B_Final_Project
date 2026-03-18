# RP1B_Final_Project

Included in this repository are the three Python scripts that were used to generate the results discussed in this project.

## rp1b_script.py

`rp1b_script.py` was initially run on the data included in the `refSeqs_and_VCF_data` directory.

The script takes the following arguments as input from the command line:
- `-i` is the input multisample VCF file path
- `-r` is the input reference genome file path (FASTA format for example)
- `-j` is the job title which is used to distinguish files generated when running on multiple occasions
- `-w` can be optionally used to alter the window size (measured in Kb) used during the identification of recombinant SNPs (the default value is 2, representing a window size of 2000 bases)
- `-s` can be optionally used to alter the threshold of SNPs required in the given window when deciding whether to classify the region as recombinant or not

## rp1b_plots.py

`rp1b_plots.py` plots the graphs shown as Figures 1-5 and Figure 9 of my written report.

This script was designed to be run specifically on the data produced during this project from `rp1b_script.py`.  Therefore, no arguments are required.  **However**, the job title used during the running of `rp1b_script.py` must be consistent with the file names used here.  **Therefore** it is suggested that when running `rp1b_script.py` on the data again, the same job titles are used and that the commands follow what is advised previously.

## rp1b_dim_reduction.py

`rp1b_dim_reduction.py` runs PCA, UMAP, and t-SNE before plotting the graphs shown as Figures 6-8 of my written report.

Similarly to `rp1b_plots.py`, this script was designed to be run specifically on the data produced when running `rp1b_script.py` on the data provided in the `refSeqs_and_VCF_data` directory.

## Running all three scripts

As was done during this project, all three scripts were run in the order shown using the commands listed below:

**Note**: all files were located within the same directory.

For *S. typhimurium*: `python rp1b_script.py -i parsnp.vcf -r S_typhimurium_HC5_296366_SAL_TC7523AA_AS.result.fasta.ref -j typh`

For *S. agalactiae*: `python rp1b_script.py -i s_agalactiae_parsnp.vcf -r GCF_015221735.2_ASM1522173v2_genomic.fna -j agalactiae`

For *S. epidermidis*: `python rp1b_script.py -i s_epidermidis_parsnp.vcf -r GCF_000011925.1_ASM1192v1_genomic.fna -j epidermidis`

For creating Figures 1-5 & 9: `python rp1b_plots.py`

For creating Figures 6-8: `python rp1b_dim_reduction.py`

## Parsnp

To obtain the multi sample VCF files used in this study, Parsnp was run on bacterial populations of 1000 genomes.  An example of the command used to generate the multi sample VCF for *S. typhimurium* is shown below:

`parsnp -r ../../team_genomes/reference/S_typhimurium_HC5_296366_SAL_TC7523AA_AS.result.fasta -d ../..team_genomes/fasta -o parsnp_output -p 8 --vcf`

## Conda Environment

The conda environment with all required packages, as well as the correct installation of Parsnp can be found at the following: `rp1b_env.yaml`
