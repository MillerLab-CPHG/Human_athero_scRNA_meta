library(Seurat)
library(SeuratDisk)
library(scRNAutils)
library(tidyverse)
library(ggsci)

###########################################################
# This script will be a pilot run for the scDRS analyses. #
###########################################################

# First we need to convert a seurat object into an. h5ad file.
files_dir = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/Wirka_scRNA_files/"
wirka_data = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Wirka_human_coronaries_scRNA_data/past_analyses/rds_objects/scRNA_athero_data_seurat_manual_labels1_res1.7.rds")
wirka_h5_seurat_file = paste(files_dir, "Wirka_data.h5Seurat", sep = "")

# Convert to .h5Seurat file
SaveH5Seurat(wirka_data, filename=wirka_h5_seurat_file)
Convert(wirka_h5_seurat_file, dest = "h5ad", assay = "RNA")


###########################################################################################
# Looks like there's an issue with negative values in normalized data during the conversion
# Authors suggest to create the .h5ad files using just RNA raw counts 

# Get sparse matrix
wirka_sparse_mt = wirka_data@assays$RNA@counts
wirka_new_seurat_obj = CreateSeuratObject(wirka_sparse_mt)
wirka_raw_h5_seurat_file = paste(files_dir, "Wirka_raw_data.h5Seurat", sep="")

# Convert to .h5Seurat and .h5ad files
SaveH5Seurat(wirka_new_seurat_obj, filename=wirka_raw_h5_seurat_file)
Convert(wirka_raw_h5_seurat_file, dest = "h5ad")