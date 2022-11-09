#!/bin/bash

# Source conda environment that has the correct package versions for GWAS munging
#conda activate munge_ldsc

# Run the mtag_munge script to make sure that GWAS summary stats are correctly formatted for LDSC
python /project/cphg-millerlab/Jose/CELLECT/ldsc/mtag_munge.py \
--sumstats /project/cphg-millerlab/Jose/GWAS_summary_stats/CAD_van_der_Harst/van_der_Harst_snATAC_manuscript/CAD_summary_for_ldsc.txt.gz \
--merge-alleles /project/cphg-millerlab/Jose/CELLECT/data/ldsc/w_hm3.snplist \
--out /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/munged_GWAS_summary_stats/CAD_GWAS_van_der_Harst/munged_CAD_van_der_Harst_snATAC_manuscript/munged_CAD_van_der_Harst \
--n-value 296525 \
--keep-pval \
--p PVAL \
--chunksize 500000
