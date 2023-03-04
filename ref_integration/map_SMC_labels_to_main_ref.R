library(tidyverse)
library(data.table)



# How many SMC barcodes do we have in the main ref (29024)
SMCs_main_ref = rpca_int_sct_v3@meta.data[rpca_int_sct_v3@meta.data$level1_annotations == "SMC",]
SMC_main_ref_barcodes = rownames(SMCs_main_ref)
head(SMC_main_ref_barcodes)
length(SMC_main_ref_barcodes)

# How many SMC barcodes do we have in the subclustered data (27605)
#SMC_subset_anno = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte", "Fibrochondrocyte",
#                    "SMC2", "SMC3", "Foam-like")
#SMC_subset = rpca_smc_fibro_subset_v3@meta.data[rpca_smc_fibro_subset_v3@meta.data$prelim_annotations %in% SMC_subset_anno, ]

SMC_subset_barcodes = rownames(rpca_smc_fibro_subset_v3@meta.data)
length(SMC_subset_barcodes)

# How many barcodes match
length(setdiff(SMC_main_ref_barcodes, SMC_subset_barcodes))

# Create a new vector where we we'll update the SMC annotations
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


# Add new annotations into metadata
rpca_int_sct_v3@meta.data$new_level2_anno = names(new_level2_barcodes)
DimPlot(rpca_int_sct_v3, group.by = "new_level2_anno", raster = FALSE, label = TRUE) & custom_theme & npg_scale2 + 
  theme(legend.position = "none")

FeaturePlot(rpca_int_sct_v3, features = c("LMOD1", "MYH11", "CNN1", "ACTA2"), raster = FALSE, order = TRUE) & custom_theme & new_scale3






