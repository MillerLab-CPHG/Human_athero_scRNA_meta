library(Seurat)
library(tidyverse)
library(data.table)
library(broom)
library(ggrepel)


# Set seed for reproducibility
set.seed(1)

# Source our own scRNA analysis utils functions 
source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")

##########################################################################################
# This script will calculate correlations for a gene of interest with all of other genes #
# for cell types of interest.                                                            #
##########################################################################################

# Load seurat object with pericytes, SMCs and Fibroblasts
rpca_smc_fibro_subset_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v1.rds")
DefaultAssay(rpca_smc_fibro_subset_v3) = "SCT"

# Get sct counts from matrix
rpca_smc_fibrosct_counts = rpca_smc_fibro_subset_v3@assays$SCT@data
class(rpca_smc_fibrosct_counts)

# Get metadata
rpca_smc_fibro_metadata = rpca_smc_fibro_subset_v3@meta.data

# Define cell types and target genes for correlations
cell_types = c("Contractile_SMC", "Fibrochondrocyte")
genes_of_interest = c("MYH11", "CNN1", "TAGLN", "LMOD1", "TPM2", "ACTA2", "PCOLCE2", "OMD",
                      "IBSP", "COMP", "LUM", "DCN", "VCAM1", "KLF4", "MGP",
                      "SERPINE2", "CYTL1", "SOX9")


# Calculate correlations for CRTAC1
# For this, we will source our custom functions (calc_gene_cors and plot_cors) 
# from the Utils script
smc_fibrochondro_crtac1 = calc_gene_cors(cell_types = cell_types,
                                         exp_matrix = rpca_smc_fibrosct_counts,
                                         metadata_df = rpca_smc_fibro_metadata, 
                                         target_gene = "CRTAC1", 
                                         cor_method = "pearson")


# Plot gene correlations for CRTAC1
crtac1_cors_plot = plot_cors(smc_fibrochondro_crtac1, "CRTAC1", "SMCs and Fibrochondrocytes", 
                             genes_of_interest) + 
  custom_theme +
  ggtitle("") + 
  theme(aspect.ratio = 1.5)
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure4/Fig4e_CRTAC1_pearson_cors.pdf",
       plot = crtac1_cors_plot, width = 6, height = 8)
  

###########################################################
# Plot correlations between FHL5 and calcification drivers

# Get barcodes for SMCs and fibrochondrocytes 
smc_barcodes = rpca_smc_fibro_metadata %>% 
  filter(prelim_annotations %in% cell_types) %>%
  rownames()

# Subset matrix to keep only SMCs and fibrochondrocytes 
smc_fibrochondro_counts = rpca_smc_fibrosct_counts[, smc_barcodes]
smc_fibrochondro_counts_df = t(as.data.frame(smc_fibrochondro_counts))
smc_fibrochondro_counts_df[1:4, 1:4]

as.data.frame(smc_fibrochondro_counts_df) %>%
  ggplot(aes(x=FHL5, y=COL1A2)) + 
  geom_smooth(method = "lm") + 
  geom_point(size=1) + 
  custom_theme



