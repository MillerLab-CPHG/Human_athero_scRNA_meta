library(Seurat)
library(SeuratDisk)
library(scRNAutils)
library(tidyverse)
library(ggsci)

#########################################################################################################
# This script will produce the files necessary for the scDRS analyses with the meta-analyzed reference. #
#########################################################################################################


# Load reference
rpca_int_sct_v3_1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_clustered_v3_1.rds")
unique(rpca_int_sct_v3_1$updated_level2_annotations)
DimPlot(rpca_int_sct_v3_1, group.by = "updated_level2_annotations", raster = FALSE) + custom_theme()

###############################################################
# First we need to convert a seurat object into an. h5ad file.
# We should have 2 .h5ad files. One with all the normalized data and one with only the raw counts
files_dir = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/whole_ref_scRNA_files/"
rpca_int_sct_h5_seurat_file = paste(files_dir, "rpca_int_sct_v3_1.h5Seurat", sep = "")

# Convert to .h5Seurat and .h5ad files. Store counts from SCT assay
SaveH5Seurat(rpca_int_sct_v3_1, filename=rpca_int_sct_h5_seurat_file)
Convert(rpca_int_sct_h5_seurat_file, dest = "h5ad", assay = "SCT")


###########################################################################################
# Looks like there's an issue with negative values in normalized data during the conversion
# Authors suggest to create the .h5ad files using just RNA raw counts 

# Get sparse matrix
rpca_int_sct_sparse_mt = rpca_int_sct_v3_1@assays$RNA@counts
rpca_int_new_seurat_obj = CreateSeuratObject(rpca_int_sct_sparse_mt)
rpca_int_raw_h5_seurat_file = paste(files_dir, "rpca_int_raw_data.h5Seurat", sep="")

# Convert to .h5Seurat and .h5ad files
SaveH5Seurat(rpca_int_new_seurat_obj, filename=rpca_int_raw_h5_seurat_file)
Convert(rpca_int_raw_h5_seurat_file, dest = "h5ad")


# Generate a covariate file to run the scDRS enrichment 
# What do we include as covariates? N genes, sex, lesion category, arterial origin
metadata = rpca_int_sct_v3_1@meta.data

new_meta_df = metadata %>%
  select(sex, arterial_origin, 
         sample_disease_status)
write.table(new_meta_df,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/categorical_covariate_file.cov",
            quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

# Need to define a scheme for converting categories into numerical values
# Non-lesion=0; Lesion=1
# Male=0; Female=1
# Carotid=0; Coronary=1

cov_file = new_meta_df %>%
  mutate(sex = case_when(sex == "males" ~ 0, 
                         TRUE ~ 1),
         sample_disease_status = case_when(sample_disease_status == "lesion" ~ 1,
                                           TRUE ~ 0),
         arterial_origin = case_when(arterial_origin == "coronary" ~ 1,
                                     TRUE ~ 0))

write.table(cov_file,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/covariate_file.cov",
            quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
  










