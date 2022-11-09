#!/bin/bash

# Run the mtag_munge script to make sure that GWAS summary stats are correctly formatted for LDSC
python /project/cphg-millerlab/Jose/CELLECT/ldsc/mtag_munge.py \
--sumstats /project/cphg-millerlab/Jose/GWAS_summary_stats/GWAS_Carotid_plaque/Plaque_summary_for_ldsc.txt.gz \
--merge-alleles /project/cphg-millerlab/Jose/CELLECT/data/ldsc/w_hm3.snplist \
--out /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/munged_GWAS_summary_stats/Carotid_plaque_GWAS/munged_Carotid_plaque \
--n-value 48000 \
--keep-pval \
--p PVAL \
--chunksize 500000
