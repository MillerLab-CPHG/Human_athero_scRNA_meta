#!/bin/bash

# Source conda environment that has the correct package versions for GWAS munging
#conda activate munge_ldsc

# Run the mtag_munge script to make sure that GWAS summary stats are correctly formatted for LDSC
python /project/cphg-millerlab/Jose/CELLECT/ldsc/mtag_munge.py \
--sumstats /project/cphg-millerlab/Jose/GWAS_summary_stats/CAD_MVP_EUR_meta/MVP_EUR_meta_CAD_NatMed2022.txt \
--merge-alleles /project/cphg-millerlab/Jose/CELLECT/data/ldsc/w_hm3.snplist \
--out /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/munged_GWAS_summary_stats/CAD_GWAS_MVP_EUR/munged_CAD_MVP \
--n-value 773268 \
--keep-pval \
--p PVAL \
--chunksize 500000
