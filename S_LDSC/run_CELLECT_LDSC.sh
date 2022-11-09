#!/bin/bash

# Load modules. Looks like python3.8 is needed for the script to finish running but need to do a few more tests. 
module load anaconda/2020.11-py3.8
module load bzip2/1.0.6
module load snakemake/6.0.5

# Run ldsc analysis using CELLECT wrapper.
# Before running this command, it might be necessary to gunzip the gene expression specificity file. Somehow it seems to help make things run more smoothly. 
snakemake --use-conda -j -s /project/cphg-millerlab/Jose/CELLECT/cellect-ldsc.snakefile --configfile cellect_config_multi_GWAS.yml
