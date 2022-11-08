#!/bin/bash

# Annotate cell types based on Immune_All_AddPIP model (no majority voting)
celltypist --indata /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/rpca_int_sct_v3_norm_counts.csv \
--model Immune_All_AddPIP.pkl \
--outdir /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/celltypist_outputs/Immune_All_AddPIP \
--prefix celltypist_All_Immune_AddPIP_ \
--transpose-input \
--plot-results
