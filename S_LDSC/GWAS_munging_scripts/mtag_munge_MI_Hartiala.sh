#!/bin/bash

# Run the mtag_munge script to make sure that GWAS summary stats are correctly formatted for LDSC
python /project/cphg-millerlab/Jose/CELLECT/ldsc/mtag_munge.py \
--sumstats /project/cphg-millerlab/Jose/GWAS_summary_stats/GWAS_MI_Hartiala/MI_Hartiala_summary_for_ldsc.txt \
--merge-alleles /project/cphg-millerlab/Jose/CELLECT/data/ldsc/w_hm3.snplist \
--out /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/munged_GWAS_summary_stats/MI_GWAS_Hartiala/munged_MI_Hartiala_meta \
--n-value 639221 \
--keep-pval \
--p PVAL \
--chunksize 500000
