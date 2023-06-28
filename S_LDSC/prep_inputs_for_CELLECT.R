library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)


################################################################################################################
# This script will prep matrix and metadata inputs for creating gene expression specificity files with CELLEX  #
################################################################################################################

######################################################################
# PREP FILES FOR WIRKA TESTS

# Prep some test files for cell expression specificity
wirka = readRDS("/project/cphg-millerlab/peak_calls_for_Gaelle/scRNA_athero_data_seurat_manual_labels1_res1.7.rds")

# Extract normalized counts matrix
norm_counts = wirka@assays$RNA@data
head(norm_counts)
write.csv(norm_counts, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/test_files/wirka_norm_counts_no_quotes.csv", 
          row.names = TRUE, quote = FALSE)

# Save metadata
cell_type_labels = wirka@meta.data$manually_annotated_labels
cell_type_metadata_df = data.frame(cell_type = cell_type_labels)
rownames(cell_type_metadata_df) = rownames(wirka@meta.data)
write.csv(cell_type_metadata_df,
          "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/test_files/wirka_cell_type_metadata.csv",
          row.names = TRUE)

# Switch gene symbols by Ensembl IDs
wirka_gene_symbols = rownames(norm_counts)
wirka_ensembl_ids_df = get_human_ensembl_ids(c("MYH11", "ACTA2", "VCAN"), human_biomart = human)




##############################################################################
# Extract normalized expression counts from rpca_int_sct ref v3 (118578 CELLS)
meta_athero_sct_norm_counts_matrix = rpca_int_sct_v3@assays$SCT@data
class(meta_athero_sct_norm_counts_matrix)
meta_athero_sct_counts_df = as.data.frame(meta_athero_sct_norm_counts_matrix)
rna_list = list(rna_counts_matrix)

# Write counts matrix to csv (this doesn't work for extremely large sparse matrices)
# write.csv(meta_athero_sct_norm_counts_matrix,
#           "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/rpca_int_sct_v3_sct_norm_counts_matrix.csv",
#           row.names = TRUE, quote = FALSE)

# Try write.matrix (good strategy for very large sparse matrices )
write(colnames(meta_athero_sct_norm_counts_matrix), 
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/rpca_int_sct_v3_sct_matrix_colnames.txt")

writeMM(meta_athero_sct_norm_counts_matrix, 
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/rpca_int_sct_v3_sct_sparse_matrix.txt")

# Matrix has 23381 genes x 118578 cells 
dim(rna_counts_matrix)
head(rownames(rna_counts_matrix))

# Write metadata for level 1 cell type annotations
cell_type_labels = rpca_int_sct_v3@meta.data$level1_annotations
cell_type_metadata_df = data.frame(cell_type = cell_type_labels)
rownames(cell_type_metadata_df) = rownames(rpca_int_sct_v3@meta.data)
write.csv(cell_type_metadata_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/rpca_int_sct_v3_level1_cell_type_annotations_metadata_quotes.csv",
          row.names = TRUE)


#####################################################################
# Extract normalized expression counts from smc_fibro_peri subset v3

# Write matrix for all pericytes, SMCs and fibroblasts into sparse matrix format
smc_peri_fibro_matrix = rpca_smc_fibro_subset_v3@assays$SCT@data
write(colnames(smc_peri_fibro_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/LDSC_Peri_SMC_Fibro/rpca_smc_fibro_peri_sct_matrix_colnames.txt")
write(rownames(smc_peri_fibro_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/LDSC_Peri_SMC_Fibro/rpca_smc_fibro_peri_sct_matrix_rownames.txt")
writeMM(smc_peri_fibro_matrix,
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/LDSC_Peri_SMC_Fibro/rpca_smc_fibro_peri_sct_sparse_matrix.txt")

# Write metadata for pericytes, SMCs and fibroblasts
smc_peri_fibro_labels = rpca_smc_fibro_subset_v3@meta.data$prelim_annotations
smc_peri_fibro_metadata_df = data.frame(cell_type = smc_peri_fibro_labels)
rownames(smc_peri_fibro_metadata_df) = rownames(rpca_smc_fibro_subset_v3@meta.data)
write.csv(smc_peri_fibro_metadata_df, 
          "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/LDSC_Peri_SMC_Fibro/rpca_smc_fibro_prelim_annotations_metadata.csv",
          row.names = TRUE, quote = FALSE)

#####################################################################
# Extract normalized expression counts of cells from lesions (level 1 annotations) v3

# Subset object to contain only cells from lesions (there are 59691 cells)
lesion_subset_v3 = subset(rpca_int_sct_v3_1, subset = sample_disease_status == "lesion")
lesion_subset_v3

# Write matrix for lesion cells into sparse matrix format
lesion_matrix = lesion_subset_v3@assays$SCT@data
write(colnames(lesion_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/lesion/input_matrix_files/rpca_lesion_matrix_colnames.txt")
write(rownames(lesion_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/lesion/input_matrix_files/rpca_lesion_matrix_rownames.txt")
writeMM(lesion_matrix,
        file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/lesion/input_matrix_files/rpca_lesion_sct_sparse_matrix.txt")

# Write metadata for pericytes, SMCs and fibroblasts
lesion_labels = lesion_subset_v3@meta.data$level1_annotations
lesion_metadata_df = data.frame(cell_type=lesion_labels)
rownames(lesion_metadata_df) = rownames(lesion_subset_v3@meta.data)
write.csv(lesion_metadata_df, 
          "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/lesion/input_matrix_files/rpca_lesion_level1_annotations_metadata.csv",
          row.names = TRUE, quote = FALSE)


################################################################################################
# Extract normalized expression counts of cells from non-lesion samples (level 1 annotations) v3

# Subset object to contain only cells from lesions (there are 58887 cells)
non_lesion_subset_v3 = subset(rpca_int_sct_v3_1, subset = sample_disease_status == "non_lesion")
non_lesion_subset_v3

# Write matrix for lesion cells into sparse matrix format
non_lesion_matrix = non_lesion_subset_v3@assays$SCT@data
write(colnames(non_lesion_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/non_lesion/input_matrix_files/rpca_non_lesion_matrix_colnames.txt")
write(rownames(non_lesion_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/non_lesion/input_matrix_files/rpca_non_lesion_matrix_rownames.txt")
writeMM(non_lesion_matrix,
        file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/non_lesion/input_matrix_files/rpca_non_lesion_sct_sparse_matrix.txt")

# Write metadata for pericytes, SMCs and fibroblasts
non_lesion_labels = non_lesion_subset_v3@meta.data$level1_annotations
non_lesion_metadata_df = data.frame(cell_type=non_lesion_labels)
rownames(non_lesion_metadata_df) = rownames(non_lesion_subset_v3@meta.data)
write.csv(non_lesion_metadata_df, 
          "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Lesion_status/non_lesion/input_matrix_files/rpca_non_lesion_level1_annotations_metadata.csv",
          row.names = TRUE, quote = FALSE)


##########################################################################
# Extract normalized expression counts of cells from Endothelial partition

# Define vector with names of endothelial subtypes
endo_vec = c("Inflammatory_EC", "Angiogenic/Vasa_vasorum_EC", 
             "Intimal_EC", "Lymphatic_EC", "EndoMT_EC")

# Subset object to contain only cells from the Endothelial compartment 
endo_subset_v3 = subset(rpca_int_sct_v3_1, subset = updated_level2_annotations %in% endo_vec)
endo_subset_v3

# Check annotations (we only kept Endo subtypes)
unique(endo_subset_v3@meta.data$updated_level2_annotations)

# Write matrix for lesion cells into sparse matrix format
endo_matrix = endo_subset_v3@assays$SCT@data
write(colnames(endo_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Endothelial/input_matrix_files/rpca_endo_matrix_colnames.txt")
write(rownames(endo_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Endothelial/input_matrix_files/rpca_endo_matrix_rownames.txt")
writeMM(endo_matrix,
        file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Endothelial/input_matrix_files/rpca_endo_sct_sparse_matrix.txt")

# Write metadata for pericytes, SMCs and fibroblasts
endo_labels = endo_subset_v3@meta.data$updated_level2_annotations
endo_metadata_df = data.frame(cell_type=endo_labels)
rownames(endo_metadata_df) = rownames(endo_subset_v3@meta.data)
write.csv(endo_metadata_df, 
          "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Endothelial/input_matrix_files/rpca_endo_level2_annotations_metadata.csv",
          row.names = TRUE, quote = FALSE)

# Check we added the proper annotations to metadata
unique(endo_metadata_df$cell_type)

###################################################################################################
# Extract normalized expression counts of cells from the Myeloid partition (level 1 annotations) v3

# Define vector with names of Myeloid subtypes
myeloid_vec = c("Foamy_Mac1", "Foamy_Mac2", "NAMPT_Neutrophils", "Phagocytosis_Mac",
                "Tissue_resident_Mac", "Monocytes", "Proliferating_myeloid",
                "Inflammatory_Mac", "Monocytes/DC", "cDC", "Mast_cell")

# Subset object to contain only cells from the Myeloid lineage (there are 24142 cells)
myeloid_subset_v3 = subset(rpca_int_sct_v3_1, subset = updated_level2_annotations %in% myeloid_vec)
myeloid_subset_v3

# Write matrix for lesion cells into sparse matrix format
myeloid_matrix = myeloid_subset_v3@assays$SCT@data
write(colnames(myeloid_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Myeloid/input_matrix_files/rpca_myeloid_matrix_colnames.txt")
write(rownames(myeloid_matrix),
      file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Myeloid/input_matrix_files/rpca_myeloid_matrix_rownames.txt")
writeMM(myeloid_matrix,
        file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Myeloid/input_matrix_files/rpca_myeloid_sct_sparse_matrix.txt")

# Write metadata for Myeloid subtypes
myeloid_labels = myeloid_subset_v3@meta.data$updated_level2_annotations
myeloid_metadata_df = data.frame(cell_type=myeloid_labels)
rownames(myeloid_metadata_df) = rownames(myeloid_subset_v3@meta.data)
write.csv(myeloid_metadata_df, 
          "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/Myeloid/input_matrix_files/rpca_myeloid_level2_annotations_metadata.csv",
          row.names = TRUE, quote = FALSE)

# Check we have the proper myeloid annotations
unique(myeloid_metadata_df$cell_type)


######################################################################
# Alternative: convert seurat object into loom object
#rpca_int_sct_v3_loom = as.loom(rpca_int_sct_v3, 
#                               filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/rpca_int_sct_v3.loom", 
#                               verbose = TRUE)

# Replace gene symbols in the matrix with Ensembl IDs
#matrix_gene_symbols = rownames(rna_counts_matrix)
#length(matrix_gene_symbols)

# Convert gene symbols to Ensembl IDs 
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#ensembl_ids_df = get_human_ensembl_ids(matrix_gene_symbols, human_biomart = human)
#dim(ensembl_ids_df)
#ensembl_ids_df[ensembl_ids_df$external_gene_name == "NFKBIA", ]


