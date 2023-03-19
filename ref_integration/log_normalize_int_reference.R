library(Seurat)
library(tidyverse)
library(data.table)


##################################################################################
# The goal of this script is to log normalize the entire reference so we can use #
# these values for tasks such as label transfer to snATAC-seq data.              #
##################################################################################


# Load main reference and set assay to RNA
DefaultAssay(rpca_int_sct_v3) = "RNA"

# Normalize the reference using the global-scaling normalization method.
# This normalizes gene expression for each cell by the total expression,
# multiplies this by a scale factor (10000) and log-transforms the result. 
# Log normalized values should be stored within rpca_int[["RNA"]]@data. 
rpca_int_sct_v3_1 = NormalizeData(rpca_int_sct_v3, 
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000) 


# Find variable features
rpca_int_sct_v3_1 = FindVariableFeatures(rpca_int_sct_v3_1, 
                                         selection.method = "vst",
                                         nfeatures = 2000)


# Scale the data
# Do this for only variable features. 
rpca_int_sct_v3_1 = ScaleData(rpca_int_sct_v3_1)

# Plot log normalized gene expression. How does it compare to
# SCTransform normalized expression? 
FeaturePlot(rpca_int_sct_v3_1, features = c("LMOD1", "MYH11", "CNN1", 
                                            "VCAN", "LTBP1", "IBSP"),
            order = TRUE, raster = FALSE) & custom_theme & new_scale3 

# Save rds object for rpca_int_sct v3.1
saveRDS(rpca_int_sct_v3_1, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_clustered_v3_1.rds")


