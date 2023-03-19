library(Seurat)
library(tidyverse)
library(data.table)



##########################################################################################
# The goal of this script is to take the cell annotations from the subclustered SMCs and #
# add them back into the main meta-analyzed scRNA reference.                             #
##########################################################################################

# How many SMC barcodes do we have in the main ref (29024)
SMCs_main_ref = rpca_int_sct_v3@meta.data[rpca_int_sct_v3@meta.data$level1_annotations == "SMC",]
SMC_main_ref_barcodes = rownames(SMCs_main_ref)
head(SMC_main_ref_barcodes)
length(SMC_main_ref_barcodes)

# How many SMC barcodes do we have in the subclustered data (27605)
#SMC_subset_anno = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte", "Fibrochondrocyte",
#                    "SMC2", "SMC3", "Foam-like")
#SMC_subset = rpca_smc_fibro_subset_v3@meta.data[rpca_smc_fibro_subset_v3@meta.data$prelim_annotations %in% SMC_subset_anno, ]

# Get barcodes from SMC subset
SMC_subset_barcodes = rownames(rpca_smc_fibro_subset_v3@meta.data)
length(SMC_subset_barcodes)

# Check for any potential discrepancies
length(setdiff(SMC_main_ref_barcodes, SMC_subset_barcodes))

# Create a new vector where we'll update the SMC annotations
new_level2_barcodes = rownames(rpca_int_sct_v3@meta.data)
names(new_level2_barcodes) = rpca_int_sct_v3@meta.data$level2_annotations
names(SMC_subset_barcodes) = rpca_smc_fibro_subset_v3$prelim_annotations
#smc_barcodes = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte", "Fibrochondrocyte", "SMC2", "SMC3", "Foam-like")
for (i in seq_len(length(new_level2_barcodes))) { 
  if (names(new_level2_barcodes[i]) %in% c("SMC", "Pericyte")) { 
    idx = match(new_level2_barcodes[i], SMC_subset_barcodes)
    names(new_level2_barcodes)[i] = names(SMC_subset_barcodes[idx])
    }
}

# Add new annotations into metadata and visualize UMAP embeddings
rpca_int_sct_v3@meta.data$updated_level2_annotations = names(new_level2_barcodes)
DimPlot(rpca_int_sct_v3, group.by = "updated_level2_annotations", raster = FALSE, label = TRUE, label.size = 5,
        repel = TRUE) & custom_theme + 
  theme(legend.position = "none")


# Check a few markers as a sanity check
FeaturePlot(rpca_int_sct_v3, features = c("LMOD1", "MYH11", "CNN1", "VCAN", "LTBP1", "IBSP"), raster = FALSE, order = TRUE) & custom_theme & new_scale3






